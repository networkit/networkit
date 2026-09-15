#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include <networkit/auxiliary/Log.hpp>

#include "SearchGraph.hpp"

namespace NetworKit {
namespace IsomorphismDetails {

namespace {

/**
 * Largest node id bound for which the bit matrix is still worth building.
 *
 * The CSR costs memory per edge; the matrix costs a bit per *ordered pair of node ids*, so it is
 * the one part of the snapshot whose size ignores how sparse the graph is. This bound is the
 * integer square root of 64 MiB expressed in bits, which is orders of magnitude beyond any real
 * pattern. It is a soft boundary: crossing it costs speed, not correctness, because hasEdge()
 * simply falls back to the CSR.
 */
constexpr count maxMatrixNodes = 23170;

/**
 * Sort every slice of a CSR ascending by neighbour, keeping the parallel label array aligned.
 *
 * Without labels this is a plain sort over the heads. With them the two arrays have to move
 * together: a label is identified by *where it sits*, so sorting the heads on their own would hand
 * every arc some other arc's label, silently and without any of it failing to compile or crash.
 * Sorting (head, label) pairs and unzipping is the least error-prone way to say that, and the
 * scratch buffer is reused across nodes so it costs one allocation for the whole graph.
 */
void sortSlices(const std::vector<index> &first, std::vector<node> &head, std::vector<index> &label,
                count z) {
    if (label.empty()) {
        // Sorting over pointers rather than iterators keeps the offsets unsigned.
        for (node u = 0; u < z; ++u)
            std::sort(head.data() + first[u], head.data() + first[u + 1]);
        return;
    }

    std::vector<std::pair<node, index>> slice;
    for (node u = 0; u < z; ++u) {
        const index begin = first[u];
        const index end = first[u + 1];

        slice.clear();
        for (index i = begin; i < end; ++i)
            slice.emplace_back(head[i], label[i]);

        // Lexicographic, so parallel arcs land next to each other in a defined order rather than
        // whichever one the scatter happened to write first.
        std::sort(slice.begin(), slice.end());

        for (index i = begin; i < end; ++i) {
            head[i] = slice[i - begin].first;
            label[i] = slice[i - begin].second;
        }
    }
}

/**
 * Drop self-loops and collapse parallel edges in a sorted CSR, in place.
 *
 * Both are meaningless to a subgraph search - the two matching semantics only ever constrain
 * pairs of distinct nodes - and both actively break it if left in: a repeated neighbour makes the
 * search enumerate the same candidate twice and report the same match twice, and either one
 * inflates the degrees that every feasibility rule prunes on.
 *
 * The slices are already sorted, so equal entries are adjacent and one linear scan suffices.
 * Compaction happens in place because the write cursor can never overtake the read cursor. The
 * label array is carried along for the same reason the sort has to carry it: it is indexed by
 * position, so moving a head without its label corrupts both.
 *
 * @return true if a collapsed run of equal heads held labels that were not all equal, so that the
 *         one arc left cannot stand for all of them. Always false without labels.
 */
bool compactSlices(std::vector<index> &first, std::vector<node> &head, std::vector<index> &label,
                   count z) {
    const bool labelled = !label.empty();
    bool lost = false;

    index write = 0;
    for (node u = 0; u < z; ++u) {
        const index begin = first[u];
        const index end = first[u + 1];
        first[u] = write; // safe: `begin` was read before this overwrites it

        node previous = none;
        for (index i = begin; i < end; ++i) {
            const node v = head[i];
            if (v == u)
                continue;
            if (v == previous) {
                // Equal heads are adjacent after the sort, so the label this one is being dropped
                // in favour of is simply the one just written.
                if (labelled && label[i] != label[write - 1])
                    lost = true;
                continue;
            }
            previous = v;
            if (labelled)
                label[write] = label[i];
            head[write] = v;
            ++write;
        }
    }

    first[z] = write;
    head.resize(write);
    if (labelled)
        label.resize(write);

    return lost;
}

} // namespace

SearchGraph::SearchGraph(const Graph &G, bool buildMatrix, const std::vector<index> &edgeLabels)
    : lostLabels(false), maxOut(0), maxIn(0), matrixStride(0), n(G.numberOfNodes()),
      z(G.upperNodeIdBound()), directed(G.isDirected()), hasMatrix(buildMatrix) {
    // Labels are indexed by edge id, so the scatter below reads edgeLabels[eid] for every arc.
    // Without ids that index does not exist and without enough entries it runs off the end, so
    // both are refused here rather than read past. SubgraphIsomorphism::setEdgeLabels() checks the
    // same two things; this is what keeps a snapshot built any other way honest.
    if (!edgeLabels.empty()) {
        if (!G.hasEdgeIds())
            throw std::runtime_error("SearchGraph: edge labels need a graph with edge ids - call "
                                     "indexEdges() on it first");
        if (edgeLabels.size() < G.upperEdgeIdBound())
            throw std::runtime_error(
                "SearchGraph: edge label vector is shorter than the graph's upperEdgeIdBound()");
    }

    // The matrix is a request, not a demand: when it will not fit, drop it and let hasEdge() use
    // the CSR, which is how every target snapshot already works. Note the bound is the *id* bound,
    // and removeNode() lowers neither it nor the ids above it, so a small pattern carved out of a
    // large graph still asks for the large graph's matrix. That case is a caller mistake with an
    // easy fix, so it warrants a warning; a genuinely large pattern does not, since there is
    // nothing to do differently.
    if (hasMatrix && z > maxMatrixNodes) {
        if (n <= z / 2) {
            WARN("SearchGraph: skipping the adjacency matrix - the node id bound is ", z,
                 " but only ", n,
                 " nodes exist. Compact the node ids first, e.g. with "
                 "GraphTools::getCompactedGraph(). Falling back to the CSR, which is correct but "
                 "slower.");
        } else {
            INFO("SearchGraph: skipping the adjacency matrix - it needs a bit per ordered pair of ",
                 z, " node ids, and is only meant for small patterns. Falling back to the CSR.");
        }
        hasMatrix = false;
    }

    buildCSR(G, edgeLabels);
    if (hasMatrix)
        buildAdjacencyMatrix();
}

void SearchGraph::buildCSR(const Graph &G, const std::vector<index> &edgeLabels) {
    const bool labelled = !edgeLabels.empty();

    // Count the out-degree of every node into outFirst[u + 1], then turn the counts into start
    // offsets with a prefix sum. A node that was removed keeps a degree of 0 and so ends up with
    // an empty slice - indistinguishable from an isolated node by adjacency alone, which is why
    // the same pass records which ids are nodes at all.
    outFirst.assign(z + 1, 0);
    nodeExists.assign(z, false);
    for (node u = 0; u < z; ++u) {
        if (G.hasNode(u)) {
            nodeExists[u] = true;
            outFirst[u + 1] = G.degreeOut(u);
        }
    }
    std::partial_sum(outFirst.begin(), outFirst.end(), outFirst.begin());

    // Place every neighbour into its node's slice. forEdges() visits an undirected edge once,
    // oriented u >= v, so the reverse orientation has to be added here - except for a self-loop,
    // which Graph stores only once and which degreeOut() therefore also counted only once.
    //
    // The four-argument overload hands out the edge id alongside the arc, which is what lets the
    // label land at the same offset as its head in the very same pass. An undirected edge writes
    // its one id into both endpoints' slices; a directed mutual pair is two ids in two different
    // slices, so per-arc labels stay well defined either way.
    outHead.resize(outFirst[z]);
    if (labelled)
        outLabel.resize(outFirst[z], none);
    std::vector<index> cursor = outFirst;
    G.forEdges([&](node u, node v, edgeweight, edgeid eid) {
        if (labelled)
            outLabel[cursor[u]] = edgeLabels[eid];
        outHead[cursor[u]++] = v;
        if (!directed && u != v) {
            if (labelled)
                outLabel[cursor[v]] = edgeLabels[eid];
            outHead[cursor[v]++] = u;
        }
    });

    // hasEdge() binary-searches a slice and the feasibility rules intersect two of them, so both
    // rely on the order.
    sortSlices(outFirst, outHead, outLabel, z);
    lostLabels |= compactSlices(outFirst, outHead, outLabel, z);

    // Only now are the slices the sets the search reasons about, so only now is the maximum the
    // number a caller may compare a pattern degree against.
    for (node u = 0; u < z; ++u) {
        maxOut = std::max(maxOut, outDegree(u));
    }

    // For an undirected graph inBegin()/inEnd() fall back to the out-arrays, so a second copy
    // would only waste memory.
    if (directed) {
        inFirst.assign(z + 1, 0);
        for (node u = 0; u < z; ++u) {
            if (G.hasNode(u)) {
                inFirst[u + 1] = G.degreeIn(u);
            }
        }
        std::partial_sum(inFirst.begin(), inFirst.end(), inFirst.begin());

        inHead.resize(inFirst[z]);
        if (labelled)
            inLabel.resize(inFirst[z], none);
        std::vector<index> inCursor = inFirst;
        G.forEdges([&](node u, node v, edgeweight, edgeid eid) {
            if (labelled)
                inLabel[inCursor[v]] = edgeLabels[eid];
            inHead[inCursor[v]++] = u;
        });

        sortSlices(inFirst, inHead, inLabel, z);
        lostLabels |= compactSlices(inFirst, inHead, inLabel, z);

        for (node u = 0; u < z; ++u) {
            maxIn = std::max(maxIn, inDegree(u));
        }
    }
}

void SearchGraph::buildAdjacencyMatrix() {
    matrixStride = tlx::div_ceil(z, 64);
    matrix.assign(static_cast<std::size_t>(z) * matrixStride, 0);

    // Filled from the compacted CSR, not from `Graph`. Parallel edges and self-loops are already
    // gone by the time this runs, so neither needs handling here - and, more to the point, neither
    // *can* be handled differently than the CSR handled it, which is what the two hasEdge()
    // backends agreeing depends on. A removed id has an empty slice and so contributes no bits.
    for (node u = 0; u < z; ++u)
        for (const node *it = outBegin(u); it != outEnd(u); ++it)
            matrix[static_cast<std::size_t>(u) * matrixStride + *it / 64] |= uint64_t{1}
                                                                             << (*it % 64);
}

} // namespace IsomorphismDetails
} // namespace NetworKit
