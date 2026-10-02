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

/// Largest node id bound for which the adjacency matrix takes at most 64 MiB.
constexpr count maxMatrixNodes = 23170;

void sortSlices(const std::vector<index> &first, std::vector<node> &head, std::vector<index> &label,
                count z) {
    if (label.empty()) {
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

        // Sorting (head, label) pairs orders parallel arcs by label.
        std::sort(slice.begin(), slice.end());

        for (index i = begin; i < end; ++i) {
            head[i] = slice[i - begin].first;
            label[i] = slice[i - begin].second;
        }
    }
}

/// Drops self-loops and collapses parallel edges in a sorted CSR. Returns true if collapsed
/// parallel arcs had different labels.
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
                // Compare with the label of the arc that was kept for this head.
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
    if (!edgeLabels.empty()) {
        if (!G.hasEdgeIds())
            throw std::runtime_error("SearchGraph: edge labels need a graph with edge ids - call "
                                     "indexEdges() on it first");
        if (edgeLabels.size() < G.upperEdgeIdBound())
            throw std::runtime_error(
                "SearchGraph: edge label vector is shorter than the graph's upperEdgeIdBound()");
    }

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

    outFirst.assign(z + 1, 0);
    nodeExists.assign(z, false);
    for (node u = 0; u < z; ++u) {
        if (G.hasNode(u)) {
            nodeExists[u] = true;
            outFirst[u + 1] = G.degreeOut(u);
        }
    }
    std::partial_sum(outFirst.begin(), outFirst.end(), outFirst.begin());

    // forEdges() visits an undirected edge once, so this adds the reverse arc, except for a
    // self-loop, which degreeOut() counts once.
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

    sortSlices(outFirst, outHead, outLabel, z);
    lostLabels |= compactSlices(outFirst, outHead, outLabel, z);

    for (node u = 0; u < z; ++u) {
        maxOut = std::max(maxOut, outDegree(u));
    }

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

    // The matrix is filled from the compacted CSR, so both hasEdge() backends agree.
    for (node u = 0; u < z; ++u)
        for (const node *it = outBegin(u); it != outEnd(u); ++it)
            matrix[static_cast<std::size_t>(u) * matrixStride + *it / 64] |= uint64_t{1}
                                                                             << (*it % 64);
}

} // namespace IsomorphismDetails
} // namespace NetworKit
