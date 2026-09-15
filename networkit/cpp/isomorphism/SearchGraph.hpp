#ifndef NETWORKIT_CPP_ISOMORPHISM_SEARCH_GRAPH_HPP_
#define NETWORKIT_CPP_ISOMORPHISM_SEARCH_GRAPH_HPP_

// Private header of the isomorphism module. Not installed, not part of the public API.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <tlx/math/div_ceil.hpp>
#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace IsomorphismDetails {

/**
 * @return the number of values that occur in both of the sorted, strictly ascending ranges
 * `[aBegin, aEnd)` and `[bBegin, bEnd)`.
 */
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
 * A read-only snapshot of a Graph, laid out for the inner loops of a subgraph search.
 *
 * `Graph::hasEdge()` scans an adjacency list, which is too slow for a search that asks it millions
 * of times. The snapshot therefore stores:
 *
 * - CSR adjacency, separately for out- and in-arcs. The out-neighbours of node @a u are
 *   `outHead[outFirst[u] .. outFirst[u + 1])`, sorted in ascending order, so @ref hasEdge() is a
 *   binary search. An undirected snapshot stores only the out-arrays.
 * - Optionally a bit-packed adjacency matrix, which answers @ref hasEdge() in constant time. It
 *   needs `upperNodeIdBound()^2` bits, so it is only built for small graphs such as the pattern.
 * - Optionally one edge label per arc, at the same offset as the arc's head.
 * - Which ids are nodes, see @ref hasNode(), and the maximum out- and in-degree.
 *
 * The snapshot is the simple graph underlying @a G: parallel edges are collapsed and self-loops
 * are dropped, so degrees count distinct neighbours. This keeps a search from trying the same
 * candidate twice and from pruning on inflated degrees. Collapsing parallel edges with different
 * labels loses a label, which @ref collapsedLabelledEdges() reports.
 *
 * The snapshot is immutable once built, so @ref ParallelRI shares it between all workers.
 */
class SearchGraph {

public:
    /**
     * Builds the snapshot.
     *
     * @param G The graph to snapshot. Not stored.
     * @param buildMatrix Whether to build the adjacency matrix. Ignored if @ref upperNodeIdBound()
     * is too large; see @ref hasAdjacencyMatrix().
     * @param edgeLabels One label per edge of @a G, indexed by edge id, or empty for no labels.
     * @throws std::runtime_error if @a edgeLabels is not empty and @a G has no edge ids or fewer
     * than `upperEdgeIdBound()` labels are given.
     */
    SearchGraph(const Graph &G, bool buildMatrix, const std::vector<index> &edgeLabels = {});

    /// Number of nodes.
    count numberOfNodes() const noexcept { return n; }

    /// One past the largest node id. All arrays are sized by this bound.
    count upperNodeIdBound() const noexcept { return z; }

    /**
     * Whether node @a u exists. A removed node and an isolated node both have empty slices, so
     * only this tells them apart. @a u must be smaller than @ref upperNodeIdBound().
     */
    bool hasNode(node u) const noexcept { return nodeExists[u]; }

    bool isDirected() const noexcept { return directed; }

    /// Largest @ref outDegree() over all nodes, or 0 if there are none.
    count maxOutDegree() const noexcept { return maxOut; }

    /// Largest @ref inDegree() over all nodes. Equals @ref maxOutDegree() for undirected graphs.
    count maxInDegree() const noexcept { return directed ? maxIn : maxOut; }

    /// Whether the adjacency matrix was built.
    bool hasAdjacencyMatrix() const noexcept { return hasMatrix; }

    /// First out-neighbour of @a u. The range up to @ref outEnd() is sorted in strictly ascending
    /// order and does not contain @a u.
    const node *outBegin(node u) const noexcept { return outHead.data() + outFirst[u]; }

    /// One past the last out-neighbour of @a u.
    const node *outEnd(node u) const noexcept { return outHead.data() + outFirst[u + 1]; }

    /// First in-neighbour of @a u. Same as @ref outBegin() for undirected graphs.
    const node *inBegin(node u) const noexcept {
        return directed ? inHead.data() + inFirst[u] : outBegin(u);
    }

    /// One past the last in-neighbour of @a u.
    const node *inEnd(node u) const noexcept {
        return directed ? inHead.data() + inFirst[u + 1] : outEnd(u);
    }

    /// Number of distinct out-neighbours of @a u other than @a u.
    count outDegree(node u) const noexcept { return outFirst[u + 1] - outFirst[u]; }

    /// Number of distinct in-neighbours of @a u other than @a u.
    count inDegree(node u) const noexcept {
        return directed ? inFirst[u + 1] - inFirst[u] : outDegree(u);
    }

    /**
     * Whether the arc @a u -> @a v exists. Constant time with the adjacency matrix, otherwise a
     * binary search. Always false for `u == v`. Both nodes must be smaller than
     * @ref upperNodeIdBound().
     */
    bool hasEdge(node u, node v) const noexcept {
        if (hasMatrix)
            return (matrix[static_cast<std::size_t>(u) * matrixStride + v / 64] >> (v % 64)) & 1u;

        return std::binary_search(outBegin(u), outEnd(u), v);
    }

    /// Number of common out-neighbours of @a u and @a v.
    count commonOutNeighbors(node u, node v) const noexcept {
        return intersectionSize(outBegin(u), outEnd(u), outBegin(v), outEnd(v));
    }

    /// Whether the snapshot has edge labels.
    bool hasEdgeLabels() const noexcept { return !outLabel.empty(); }

    /**
     * Whether collapsing parallel edges discarded a label, that is, whether parallel edges had
     * different labels. An algorithm that supports edge labels must refuse such input. Always
     * false without labels.
     */
    bool collapsedLabelledEdges() const noexcept { return lostLabels; }

    /**
     * The label of the arc @a u -> @a v, or @ref none if the arc does not exist or the snapshot
     * has no labels. A binary search over the out-neighbours of @a u, also with the adjacency
     * matrix. The two arcs of a directed mutual pair have independent labels.
     */
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

    /// Labels of the out-arcs of @a u, in the order of the @ref outBegin() range, or nullptr
    /// without labels.
    const index *outLabelBegin(node u) const noexcept {
        return outLabel.empty() ? nullptr : outLabel.data() + outFirst[u];
    }

    /// Labels of the in-arcs of @a u, in the order of the @ref inBegin() range, or nullptr without
    /// labels. Same as @ref outLabelBegin() for undirected graphs.
    const index *inLabelBegin(node u) const noexcept {
        if (!directed)
            return outLabelBegin(u);
        return inLabel.empty() ? nullptr : inLabel.data() + inFirst[u];
    }

private:
    /**
     * Fills the CSR arrays, the edge labels, `nodeExists`, `maxOut` and `maxIn`. Scatters the arcs
     * of @a G into their slices, sorts every slice together with its labels, and then drops
     * self-loops and collapses parallel edges.
     *
     * @param G The graph to snapshot.
     * @param edgeLabels Labels indexed by edge id, or empty for no labels.
     */
    void buildCSR(const Graph &G, const std::vector<index> &edgeLabels);

    /**
     * Fills the adjacency matrix from the CSR arrays, so it must run after @ref buildCSR(). Row
     * `u` has `matrixStride` 64-bit words, and bit `v` is set iff the arc `u -> v` exists.
     */
    void buildAdjacencyMatrix();

    /// CSR out-arcs: outHead[outFirst[u] .. outFirst[u + 1]) are the out-neighbours of u.
    std::vector<index> outFirst;
    std::vector<node> outHead;

    /// CSR in-arcs. Empty for undirected graphs.
    std::vector<index> inFirst;
    std::vector<node> inHead;

    /// One label per arc, at the same offset as its head. Empty without labels; inLabel is also
    /// empty for undirected graphs.
    std::vector<index> outLabel;
    std::vector<index> inLabel;

    /// Whether collapsing parallel arcs discarded a label.
    bool lostLabels;

    /// Whether each id is a node.
    std::vector<bool> nodeExists;

    /// Maximum out- and in-degree. maxIn is 0 for undirected graphs.
    count maxOut;
    count maxIn;

    /// Bit-packed adjacency matrix, row-major with matrixStride words per row. Empty if not built.
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
