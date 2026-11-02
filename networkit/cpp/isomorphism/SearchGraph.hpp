#ifndef NETWORKIT_CPP_ISOMORPHISM_SEARCH_GRAPH_HPP_
#define NETWORKIT_CPP_ISOMORPHISM_SEARCH_GRAPH_HPP_

#include <algorithm>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit::IsomorphismDetails {

/// @ref none is a wildcard on either side.
inline bool labelsCompatible(index patternLabel, index targetLabel) noexcept {
    return patternLabel == none || targetLabel == none || patternLabel == targetLabel;
}

/// Empty label vectors put no constraint on the nodes.
inline bool nodeLabelsCompatible(std::span<const index> patternNodeLabels,
                                 std::span<const index> targetNodeLabels, node pu, node tv) {
    return patternNodeLabels.empty()
           || labelsCompatible(patternNodeLabels[pu], targetNodeLabels[tv]);
}

/**
 * A read-only snapshot of a Graph for the inner loops of a subgraph search, since
 * `Graph::hasEdge()` is too slow there. It stores sorted CSR adjacency for out- and in-arcs,
 * optionally one edge label per arc, and optionally a bit-packed adjacency matrix for small graphs
 * such as the pattern. The snapshot collapses parallel edges and drops self-loops, so degrees count
 * distinct neighbors.
 */
class SearchGraph {

public:
    /// @a edgeLabels is indexed by edge id, or empty for no labels. The snapshot skips the
    /// adjacency matrix if the node id bound of @a G is too large.
    SearchGraph(const Graph &G, bool buildMatrix, std::span<const index> edgeLabels = {});

    count numberOfNodes() const noexcept { return n; }

    count upperNodeIdBound() const noexcept { return z; }

    /// A removed node and an isolated node both have empty slices; only this tells them apart.
    bool hasNode(node u) const noexcept { return nodeExists[u]; }

    bool isDirected() const noexcept { return directed; }

    count maxOutDegree() const noexcept { return maxOut; }

    count maxInDegree() const noexcept { return directed ? maxIn : maxOut; }

    bool hasAdjacencyMatrix() const noexcept { return hasMatrix; }

    count outDegree(node u) const noexcept { return outFirst[u + 1] - outFirst[u]; }

    count inDegree(node u) const noexcept {
        return directed ? inFirst[u + 1] - inFirst[u] : outDegree(u);
    }

    /// The distinct out-neighbors of @a u other than @a u, in ascending order.
    std::span<const node> outNeighbors(node u) const noexcept {
        return {outHead.data() + outFirst[u], outDegree(u)};
    }

    /// Like @ref outNeighbors() for the in-neighbors. Equal to it for an undirected graph.
    std::span<const node> inNeighbors(node u) const noexcept {
        if (!directed)
            return outNeighbors(u);
        return {inHead.data() + inFirst[u], inDegree(u)};
    }

    bool hasEdge(node u, node v) const noexcept {
        if (hasMatrix)
            return (matrix[matrixWord(u, v)] >> (v % 64)) & 1u;

        return std::ranges::binary_search(outNeighbors(u), v);
    }

    bool hasEdgeLabels() const noexcept { return !outLabel.empty(); }

    /// Whether collapsing parallel edges with different labels discarded a label. The algorithms
    /// throw from run() for such input.
    bool collapsedLabelledEdges() const noexcept { return lostLabels; }

    /// The label of the arc @a u -> @a v, or @ref none if the arc does not exist or the snapshot
    /// has no labels.
    index edgeLabel(node u, node v) const noexcept {
        if (outLabel.empty())
            return none;

        const std::span<const node> neighbors = outNeighbors(u);
        const auto found = std::ranges::lower_bound(neighbors, v);
        if (found == neighbors.end() || *found != v)
            return none;

        return outLabels(u)[static_cast<index>(found - neighbors.begin())];
    }

    /// The labels of the arcs to @ref outNeighbors(), position for position, or an empty span
    /// without labels.
    std::span<const index> outLabels(node u) const noexcept {
        if (outLabel.empty())
            return {};
        return {outLabel.data() + outFirst[u], outDegree(u)};
    }

    /// Like @ref outLabels() for the arcs from @ref inNeighbors().
    std::span<const index> inLabels(node u) const noexcept {
        if (!directed)
            return outLabels(u);
        if (inLabel.empty())
            return {};
        return {inLabel.data() + inFirst[u], inDegree(u)};
    }

private:
    void buildCSR(const Graph &G, std::span<const index> edgeLabels);

    void buildAdjacencyMatrix();

    /// The word of the adjacency matrix that holds the bit of the arc @a u -> @a v.
    index matrixWord(node u, node v) const noexcept { return u * matrixStride + v / 64; }

    std::vector<index> outFirst;
    std::vector<node> outHead;

    /// Empty for undirected graphs, which use the out-arrays.
    std::vector<index> inFirst;
    std::vector<node> inHead;

    std::vector<index> outLabel;
    std::vector<index> inLabel;

    bool lostLabels = false;

    std::vector<bool> nodeExists;

    count maxOut = 0;
    count maxIn = 0;

    /// Row-major, with matrixStride 64-bit words per row. Empty if not built.
    std::vector<uint64_t> matrix;
    count matrixStride = 0;

    count n;
    count z;
    bool directed;
    bool hasMatrix;
};

/// The checks of a candidate that do not depend on the partial mapping: @a tv needs at least the
/// out- and in-degree of @a pu, and a compatible node label.
inline bool degreesAndLabelsFit(const SearchGraph &pattern, const SearchGraph &target,
                                std::span<const index> patternNodeLabels,
                                std::span<const index> targetNodeLabels, node pu, node tv) {
    if (target.outDegree(tv) < pattern.outDegree(pu))
        return false;

    // The in-degrees of an undirected graph are its out-degrees.
    if (pattern.isDirected() && target.inDegree(tv) < pattern.inDegree(pu))
        return false;

    return nodeLabelsCompatible(patternNodeLabels, targetNodeLabels, pu, tv);
}

/// Throws `std::runtime_error` if @a pattern or @a target lost an edge label when it collapsed
/// parallel edges. The message starts with @a algorithmName.
void refuseLostLabels(const SearchGraph &pattern, const SearchGraph &target,
                      std::string_view algorithmName);

} // namespace NetworKit::IsomorphismDetails

#endif // NETWORKIT_CPP_ISOMORPHISM_SEARCH_GRAPH_HPP_
