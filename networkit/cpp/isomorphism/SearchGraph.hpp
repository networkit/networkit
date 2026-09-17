#ifndef NETWORKIT_CPP_ISOMORPHISM_SEARCH_GRAPH_HPP_
#define NETWORKIT_CPP_ISOMORPHISM_SEARCH_GRAPH_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <tlx/math/div_ceil.hpp>
#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace IsomorphismDetails {

/// Number of values that occur in both sorted, strictly ascending ranges.
inline count intersectionSize(const node *aBegin, const node *aEnd, const node *bBegin,
                              const node *bEnd) noexcept {
    count common = 0;
    while (aBegin != aEnd && bBegin != bEnd) {
        if (*aBegin < *bBegin) {
            ++aBegin;
        } else if (*bBegin < *aBegin) {
            ++bBegin;
        } else {
            ++common;
            ++aBegin;
            ++bBegin;
        }
    }
    return common;
}

/**
 * A read-only snapshot of a Graph for the inner loops of a subgraph search, since
 * `Graph::hasEdge()` is too slow there. It stores sorted CSR adjacency for out- and in-arcs,
 * optionally one edge label per arc, and optionally a bit-packed adjacency matrix for small graphs
 * such as the pattern. The snapshot collapses parallel edges and drops self-loops, so degrees count
 * distinct neighbours.
 */
class SearchGraph {

public:
    /// @a edgeLabels is indexed by edge id, or empty for no labels. The snapshot skips the
    /// adjacency matrix if the node id bound of @a G is too large.
    SearchGraph(const Graph &G, bool buildMatrix, const std::vector<index> &edgeLabels = {});

    count numberOfNodes() const noexcept { return n; }

    count upperNodeIdBound() const noexcept { return z; }

    /// A removed node and an isolated node both have empty slices; only this tells them apart.
    bool hasNode(node u) const noexcept { return nodeExists[u]; }

    bool isDirected() const noexcept { return directed; }

    count maxOutDegree() const noexcept { return maxOut; }

    count maxInDegree() const noexcept { return directed ? maxIn : maxOut; }

    bool hasAdjacencyMatrix() const noexcept { return hasMatrix; }

    const node *outBegin(node u) const noexcept { return outHead.data() + outFirst[u]; }

    const node *outEnd(node u) const noexcept { return outHead.data() + outFirst[u + 1]; }

    const node *inBegin(node u) const noexcept {
        return directed ? inHead.data() + inFirst[u] : outBegin(u);
    }

    const node *inEnd(node u) const noexcept {
        return directed ? inHead.data() + inFirst[u + 1] : outEnd(u);
    }

    count outDegree(node u) const noexcept { return outFirst[u + 1] - outFirst[u]; }

    count inDegree(node u) const noexcept {
        return directed ? inFirst[u + 1] - inFirst[u] : outDegree(u);
    }

    bool hasEdge(node u, node v) const noexcept {
        if (hasMatrix)
            return (matrix[static_cast<std::size_t>(u) * matrixStride + v / 64] >> (v % 64)) & 1u;

        return std::binary_search(outBegin(u), outEnd(u), v);
    }

    count commonOutNeighbors(node u, node v) const noexcept {
        return intersectionSize(outBegin(u), outEnd(u), outBegin(v), outEnd(v));
    }

    bool hasEdgeLabels() const noexcept { return !outLabel.empty(); }

    /// Whether collapsing parallel edges with different labels discarded a label. An algorithm
    /// that supports edge labels must refuse such input.
    bool collapsedLabelledEdges() const noexcept { return lostLabels; }

    /// The label of the arc @a u -> @a v, or @ref none if the arc does not exist or the snapshot
    /// has no labels.
    index edgeLabel(node u, node v) const noexcept {
        if (outLabel.empty())
            return none;

        const node *begin = outBegin(u);
        const node *end = outEnd(u);
        const node *found = std::lower_bound(begin, end, v);
        if (found == end || *found != v)
            return none;

        return outLabel[outFirst[u] + static_cast<index>(found - begin)];
    }

    /// Labels aligned with the @ref outBegin() range, or nullptr without labels.
    const index *outLabelBegin(node u) const noexcept {
        return outLabel.empty() ? nullptr : outLabel.data() + outFirst[u];
    }

    const index *inLabelBegin(node u) const noexcept {
        if (!directed)
            return outLabelBegin(u);
        return inLabel.empty() ? nullptr : inLabel.data() + inFirst[u];
    }

private:
    void buildCSR(const Graph &G, const std::vector<index> &edgeLabels);

    void buildAdjacencyMatrix();

    std::vector<index> outFirst;
    std::vector<node> outHead;

    /// Empty for undirected graphs, which use the out-arrays.
    std::vector<index> inFirst;
    std::vector<node> inHead;

    std::vector<index> outLabel;
    std::vector<index> inLabel;

    bool lostLabels;

    std::vector<bool> nodeExists;

    count maxOut;
    count maxIn;

    /// Row-major, with matrixStride 64-bit words per row. Empty if not built.
    std::vector<uint64_t> matrix;
    count matrixStride;

    count n;
    count z;
    bool directed;
    bool hasMatrix;
};

} // namespace IsomorphismDetails
} // namespace NetworKit

#endif // NETWORKIT_CPP_ISOMORPHISM_SEARCH_GRAPH_HPP_
