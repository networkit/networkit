#include <algorithm>
#include <cassert>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <tlx/math/div_ceil.hpp>
#include <networkit/auxiliary/Log.hpp>

#include "SearchGraph.hpp"

namespace NetworKit {
namespace IsomorphismDetails {

namespace {

/// Whether the adjacency matrix for the node id bound @a z fits in 64 MiB. Each row takes whole
/// 64-bit words. We divide rather than multiply, so that a huge @a z cannot overflow.
bool matrixFits(count z) {
    constexpr count maxWords = (count{64} << 20) / sizeof(uint64_t);
    return z == 0 || tlx::div_ceil(z, 64) <= maxWords / z;
}

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
    // SubgraphIsomorphism::prepareRun() has validated the labels.
    assert(edgeLabels.empty() || (G.hasEdgeIds() && edgeLabels.size() >= G.upperEdgeIdBound()));

    if (hasMatrix && !matrixFits(z)) {
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

    // forEdgesOf() visits every edge at u, so an undirected edge lands in both slices.
    const count arcs = directed ? G.numberOfEdges() : 2 * G.numberOfEdges();
    outHead.reserve(arcs);
    if (labelled)
        outLabel.reserve(arcs);

    nodeExists.assign(z, false);
    outFirst.assign(z + 1, 0);
    for (node u = 0; u < z; ++u) {
        outFirst[u] = outHead.size();
        if (!G.hasNode(u))
            continue;
        nodeExists[u] = true;
        G.forEdgesOf(u, [&](node, node v, edgeweight, edgeid eid) {
            outHead.push_back(v);
            if (labelled)
                outLabel.push_back(edgeLabels[eid]);
        });
    }
    outFirst[z] = outHead.size();

    sortSlices(outFirst, outHead, outLabel, z);
    lostLabels = compactSlices(outFirst, outHead, outLabel, z);

    for (node u = 0; u < z; ++u)
        maxOut = std::max(maxOut, outDegree(u));

    if (!directed)
        return;

    // Transposing the compacted out-arcs in tail order yields sorted in-slices. They have no
    // self-loops or parallel arcs, and they carry the labels the out-arcs kept.
    inFirst.assign(z + 1, 0);
    for (const node v : outHead)
        ++inFirst[v + 1];
    std::partial_sum(inFirst.begin(), inFirst.end(), inFirst.begin());

    inHead.resize(outHead.size());
    inLabel.resize(outLabel.size());
    std::vector<index> cursor = inFirst;
    for (node u = 0; u < z; ++u) {
        for (index i = outFirst[u]; i < outFirst[u + 1]; ++i) {
            const index slot = cursor[outHead[i]]++;
            inHead[slot] = u;
            if (labelled)
                inLabel[slot] = outLabel[i];
        }
    }

    for (node u = 0; u < z; ++u)
        maxIn = std::max(maxIn, inDegree(u));
}

void SearchGraph::buildAdjacencyMatrix() {
    matrixStride = tlx::div_ceil(z, 64);
    matrix.assign(static_cast<std::size_t>(z) * matrixStride, 0);

    // We fill the matrix from the compacted CSR, so both hasEdge() backends agree.
    for (node u = 0; u < z; ++u)
        for (const node *it = outBegin(u); it != outEnd(u); ++it)
            matrix[static_cast<std::size_t>(u) * matrixStride + *it / 64] |= uint64_t{1}
                                                                             << (*it % 64);
}

void refuseLostLabels(const SearchGraph &pattern, const SearchGraph &target,
                      const std::string &algorithmName) {
    if (!pattern.collapsedLabelledEdges() && !target.collapsedLabelledEdges())
        return;

    const std::string side = pattern.collapsedLabelledEdges() ? "pattern" : "target";
    throw std::runtime_error(algorithmName
                             + " does not support parallel edges whose edge labels disagree, as in "
                               "the "
                             + side + " - see SubgraphIsomorphism::setEdgeLabels()");
}

} // namespace IsomorphismDetails
} // namespace NetworKit
