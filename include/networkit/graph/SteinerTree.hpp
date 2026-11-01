#ifndef NETWORKIT_GRAPH_STEINER_TREE_HPP_
#define NETWORKIT_GRAPH_STEINER_TREE_HPP_

#include <compare>
#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/AdjListGraph.hpp>

namespace NetworKit {

/**
 * @ingroup graph
 * Computes a light tree that connects a given set of terminal nodes, possibly through other nodes.
 * Finding the lightest one is NP-hard. This is the approximation algorithm from "A faster
 * approximation algorithm for the Steiner problem in graphs" by K. Mehlhorn (1988): the weight of
 * the result is at most 2(1 - 1/l) times the optimal one, where l is the number of leaves of an
 * optimal tree, which is at most the number of terminals. The resulting tree is then replaced by a
 * minimum spanning tree of the subgraph induced by its nodes, and leaves that are not terminals are
 * removed. This keeps the guarantee and often makes the tree lighter. It runs in O(m log n) time.
 *
 * Unweighted graphs are treated as having unit weights. Path lengths and the total weight are
 * computed in a 64-bit integer type, or in double for floating-point weights, so that they do not
 * overflow narrow edge weight types.
 *
 * @tparam GraphT The type of the input graph.
 * @tparam TreeT The type of the resulting tree, a mutable graph type.
 */
template <class GraphT,
          class TreeT = AdjListGraph<typename GraphT::NodeT, typename GraphT::EdgeWeightT>>
class SteinerTree final : public Algorithm {

public:
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;
    // The type of sums of edge weights, such as path lengths and the total weight of the tree.
    using WeightSumT =
        std::conditional_t<std::is_floating_point_v<EdgeWeightT>,
                           std::common_type_t<EdgeWeightT, double>,
                           std::conditional_t<std::is_signed_v<EdgeWeightT>, int64_t, uint64_t>>;

    static_assert(std::is_same_v<typename TreeT::NodeT, NodeT>,
                  "SteinerTree: the tree must have the same node type as the graph.");
    static_assert(std::is_same_v<typename TreeT::EdgeWeightT, EdgeWeightT>,
                  "SteinerTree: the tree must have the same edge weight type as the graph.");

    /**
     * Constructs the algorithm for the Graph @a G.
     *
     * @param[in] G An undirected graph with non-negative edge weights.
     * @param[in] terminals The nodes to connect. Duplicates are ignored.
     * @throws std::runtime_error if @a G is directed or has negative edge weights, or if one of
     * the terminals is not in @a G.
     */
    SteinerTree(const GraphT &G, std::span<const NodeT> terminals);

    /**
     * Computes the tree for the graph and terminals passed in the constructor.
     */
    void run() override;

    /**
     * Returns whether all terminals are in the same connected component, which is always the case
     * for fewer than two terminals.
     *
     * @return True if a tree connecting all terminals exists.
     */
    bool hasTree() const {
        assureFinished();
        return connected;
    }

    /**
     * Returns the tree, a weighted graph with the same node ids as the input graph. Nodes deleted
     * from the input graph are deleted from the tree too, and nodes outside the tree have no edges.
     * Its edges have unit weight if the input graph is unweighted. It has no edges if the terminals
     * are not connected.
     *
     * @return The tree connecting the terminals.
     */
    const TreeT &getTree() const {
        assureFinished();
        return tree;
    }

    /**
     * Returns the total weight of the tree.
     *
     * @return The total weight of the edges of the tree.
     */
    WeightSumT getTotalWeight() const {
        assureFinished();
        return totalWeight;
    }

private:
    // Edges are sorted by weight, ties are broken by their endpoints. The weight is a WeightSumT
    // because bridges are weighted by path lengths.
    struct Edge {
        WeightSumT weight;
        NodeT u;
        NodeT v;

        auto operator<=>(const Edge &) const = default;
    };

    // For each node: the distance to its closest terminal, that terminal (NullNodeId if no
    // terminal is reachable), and the previous node on a shortest path from it.
    struct ClosestTerminals {
        std::vector<WeightSumT> distance;
        std::vector<NodeT> terminal;
        std::vector<NodeT> parent;
    };

    const GraphT *G;
    std::vector<NodeT> terminals;
    bool connected = false;
    WeightSumT totalWeight = 0;
    TreeT tree;

    ClosestTerminals computeClosestTerminals() const;

    // Returns the nodes of the tree found by Mehlhorn's algorithm, sets connected.
    std::vector<NodeT> mehlhornTreeNodes(const ClosestTerminals &closest);

    void addInducedMinimumSpanningTree(std::span<const NodeT> treeNodes);

    void pruneNonTerminalLeaves(std::span<const NodeT> treeNodes);
};

} /* namespace NetworKit */

#include <networkit/graph/SteinerTreeImpl.hpp>

#endif // NETWORKIT_GRAPH_STEINER_TREE_HPP_
