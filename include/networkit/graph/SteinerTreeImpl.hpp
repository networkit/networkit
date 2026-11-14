#ifndef NETWORKIT_GRAPH_STEINER_TREE_IMPL_HPP_
#define NETWORKIT_GRAPH_STEINER_TREE_IMPL_HPP_

#include <algorithm>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/VectorComparator.hpp>
#include <networkit/structures/GenericUnionFind.hpp>

namespace NetworKit {

template <class GraphT, class TreeT>
SteinerTree<GraphT, TreeT>::SteinerTree(const GraphT &G, std::span<const NodeT> terminals)
    : G(&G), terminals(terminals.begin(), terminals.end()) {
    if (G.isDirected()) {
        throw std::runtime_error("SteinerTree: the graph must be undirected.");
    }
    if constexpr (std::is_signed_v<EdgeWeightT>) {
        G.forEdges([](NodeT, NodeT, EdgeWeightT weight, edgeid) -> void {
            if (weight < 0) {
                throw std::runtime_error("SteinerTree: edge weights must not be negative.");
            }
        });
    }
    for (const NodeT terminal : this->terminals) {
        if (!G.hasNode(terminal)) {
            throw std::runtime_error("SteinerTree: the terminals must be in the graph.");
        }
    }
    std::ranges::sort(this->terminals);
    // Not std::ranges::unique, whose subrange result clang 14 cannot compile with libstdc++ 12.
    this->terminals.erase(std::unique(this->terminals.begin(), this->terminals.end()),
                          this->terminals.end());
}

template <class GraphT, class TreeT>
void SteinerTree<GraphT, TreeT>::run() {
    Aux::SignalHandler handler;
    const NodeT z = G->upperNodeIdBound();

    connected = false;
    totalWeight = 0;
    tree = TreeT(z, true);
    for (NodeT u = 0; u < z; ++u) {
        if (!G->hasNode(u)) {
            tree.removeNode(u);
        }
    }

    const ClosestTerminals closest = computeClosestTerminals();
    handler.assureRunning();

    const std::vector<NodeT> treeNodes = mehlhornTreeNodes(closest);
    if (connected) {
        handler.assureRunning();
        addInducedMinimumSpanningTree(treeNodes);
        pruneNonTerminalLeaves(treeNodes);
        // Not tree.totalEdgeWeight(), which sums edge ids if EdgeWeightT is edgeid.
        tree.forEdges([&](NodeT, NodeT, EdgeWeightT weight, edgeid) -> void {
            totalWeight += static_cast<WeightSumT>(weight);
        });
    }

    hasRun = true;
}

template <class GraphT, class TreeT>
auto SteinerTree<GraphT, TreeT>::computeClosestTerminals() const -> ClosestTerminals {
    const auto z = static_cast<std::size_t>(G->upperNodeIdBound());
    ClosestTerminals closest{std::vector<WeightSumT>(z, 0),
                             std::vector<NodeT>(z, NullNodeId<NodeT>),
                             std::vector<NodeT>(z, NullNodeId<NodeT>)};
    std::vector<WeightSumT> &distance = closest.distance;
    std::vector<NodeT> &closestTerminal = closest.terminal;
    std::vector<NodeT> &parent = closest.parent;

    // A Dijkstra search from all terminals. The heap only accepts unsigned keys. Unreached nodes
    // are recognized by their closest terminal, so no distance is reserved for them.
    using HeapKeyT = std::make_unsigned_t<NodeT>;
    tlx::d_ary_addressable_int_heap<HeapKeyT, 2, Aux::LessInVector<WeightSumT>> heap{
        Aux::LessInVector<WeightSumT>{distance}};
    heap.reserve(z);
    for (const NodeT terminal : terminals) {
        closestTerminal[terminal] = terminal;
        heap.push(static_cast<HeapKeyT>(terminal));
    }

    while (!heap.empty()) {
        const auto u = static_cast<NodeT>(heap.extract_top());
        // Handles that take the weight also take the edge id: a handle (v, weight) would get two
        // nodes if EdgeWeightT is NodeT, and (u, v, weight) would get the edge id if it is edgeid.
        G->forNeighborsOf(u, [&](NodeT, NodeT v, EdgeWeightT weight, edgeid) -> void {
            const WeightSumT newDistance = distance[u] + static_cast<WeightSumT>(weight);
            if (closestTerminal[v] == NullNodeId<NodeT> || newDistance < distance[v]) {
                distance[v] = newDistance;
                closestTerminal[v] = closestTerminal[u];
                parent[v] = u;
                heap.update(static_cast<HeapKeyT>(v));
            }
        });
    }
    return closest;
}

template <class GraphT, class TreeT>
auto SteinerTree<GraphT, TreeT>::mehlhornTreeNodes(const ClosestTerminals &closest)
    -> std::vector<NodeT> {
    const NodeT z = G->upperNodeIdBound();
    const std::vector<WeightSumT> &distance = closest.distance;
    const std::vector<NodeT> &closestTerminal = closest.terminal;
    const std::vector<NodeT> &parent = closest.parent;

    // An edge between the regions of two terminals, weighted by the length of the shortest path
    // between them through it.
    std::vector<Edge> bridges;
    G->forEdges([&](NodeT u, NodeT v, EdgeWeightT weight, edgeid) -> void {
        if (closestTerminal[u] != NullNodeId<NodeT> && closestTerminal[v] != NullNodeId<NodeT>
            && closestTerminal[u] != closestTerminal[v]) {
            bridges.push_back({distance[u] + static_cast<WeightSumT>(weight) + distance[v], u, v});
        }
    });
    std::ranges::sort(bridges);

    // The bridges of a minimum spanning tree of the regions, plus the paths from their endpoints
    // to the terminals.
    std::vector<NodeT> treeNodes;
    std::vector<Visit> inTree(z, Visit::Unvisited);
    const auto addPathToTerminal = [&](NodeT u) -> void {
        for (; u != NullNodeId<NodeT> && inTree[u] == Visit::Unvisited; u = parent[u]) {
            inTree[u] = Visit::Visited;
            treeNodes.push_back(u);
        }
    };

    GenericUnionFind<NodeT> regions(z);
    count numberOfComponents = terminals.size();
    for (const Edge &bridge : bridges) {
        if (numberOfComponents <= 1) {
            break;
        }
        const NodeT regionU = regions.find(closestTerminal[bridge.u]);
        const NodeT regionV = regions.find(closestTerminal[bridge.v]);
        if (regionU != regionV) {
            regions.merge(regionU, regionV);
            --numberOfComponents;
            addPathToTerminal(bridge.u);
            addPathToTerminal(bridge.v);
        }
    }

    connected = numberOfComponents <= 1;
    return treeNodes;
}

template <class GraphT, class TreeT>
void SteinerTree<GraphT, TreeT>::addInducedMinimumSpanningTree(std::span<const NodeT> treeNodes) {
    const NodeT z = G->upperNodeIdBound();
    std::vector<Visit> inTree(z, Visit::Unvisited);
    for (const NodeT u : treeNodes) {
        inTree[u] = Visit::Visited;
    }

    // Each edge is seen from both endpoints, and self-loops never belong to a tree.
    std::vector<Edge> edges;
    for (const NodeT u : treeNodes) {
        G->forNeighborsOf(u, [&](NodeT, NodeT v, EdgeWeightT weight, edgeid) -> void {
            if (u < v && inTree[v] == Visit::Visited) {
                edges.push_back({weight, u, v});
            }
        });
    }
    std::ranges::sort(edges);

    GenericUnionFind<NodeT> components(z);
    for (const Edge &edge : edges) {
        if (components.find(edge.u) != components.find(edge.v)) {
            components.merge(edge.u, edge.v);
            tree.addEdge(edge.u, edge.v, static_cast<EdgeWeightT>(edge.weight));
        }
    }
}

template <class GraphT, class TreeT>
void SteinerTree<GraphT, TreeT>::pruneNonTerminalLeaves(std::span<const NodeT> treeNodes) {
    std::vector<Visit> isTerminal(G->upperNodeIdBound(), Visit::Unvisited);
    for (const NodeT terminal : terminals) {
        isTerminal[terminal] = Visit::Visited;
    }
    const auto isNonTerminalLeaf = [&](NodeT u) -> bool {
        return tree.degree(u) == 1 && isTerminal[u] == Visit::Unvisited;
    };

    std::vector<NodeT> leaves;
    for (const NodeT u : treeNodes) {
        if (isNonTerminalLeaf(u)) {
            leaves.push_back(u);
        }
    }

    // The tree contains at least two terminals, so removing a leaf never leaves an isolated edge.
    while (!leaves.empty()) {
        const NodeT leaf = leaves.back();
        leaves.pop_back();
        const NodeT neighbor = tree.getIthNeighbor(leaf, 0);
        tree.removeEdge(leaf, neighbor);
        if (isNonTerminalLeaf(neighbor)) {
            leaves.push_back(neighbor);
        }
    }
}

} /* namespace NetworKit */

#endif // NETWORKIT_GRAPH_STEINER_TREE_IMPL_HPP_
