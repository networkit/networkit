#ifndef NETWORKIT_DYN_COMPONENTS_CONNECTED_COMPONENTS_IMPL_HPP_
#define NETWORKIT_DYN_COMPONENTS_CONNECTED_COMPONENTS_IMPL_HPP_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {
namespace DynConnectedComponentsDetails {

/**
 * @ingroup components
 * Determines and updates the (weakly) connected components of an undirected graph.
 */
template <bool WeaklyCC = false>
class DynConnectedComponentsImpl final : public Algorithm, public DynAlgorithm {
public:
    /**
     * Create ConnectedComponents class for Graph @a G.
     *
     * @param G The graph.
     */
    DynConnectedComponentsImpl(const Graph &G, Partition &components);

    /**
     * Finds the (weakly) connected components of the input graph.
     */
    void run() override;

    /**
     * Updates the (weakly) connected components after an edge insertion or deletion.
     *
     * @param[in]	event	The event that happened (edge deletion or
     * insertion).
     */
    void update(GraphEvent e) override;

    /**
     * Updates the (weakly) connected components after a batch of edge insertions or
     * deletions.
     *
     * @param[in] batch	A vector that contains a batch of edge insertions or
     *					deletions.
     */
    void updateBatch(const std::vector<GraphEvent> &batch) override;

private:
    const Graph *G;
    Partition *componentPtr;

    std::vector<bool> isTree;
    std::unordered_map<Edge, edgeid> edgesMap;
    std::vector<count> tmpDistances;

    void addEdge(node u, node v);
    void removeEdge(node u, node v);
    void updateTree(node u, node v);
    void indexEdges();
    void init();
};

template <bool WeaklyCC>
DynConnectedComponentsImpl<WeaklyCC>::DynConnectedComponentsImpl(const Graph &G,
                                                                 Partition &components)
    : G(&G), componentPtr(&components) {

    if (!WeaklyCC && G.isDirected())
        throw std::runtime_error("Error, connected components of directed graphs cannot be "
                                 "computed, use DynWeaklyConnectedComponents instead.");

    if (WeaklyCC && !G.isDirected())
        throw std::runtime_error("Error, weakly connected components of an undirected graph cannot "
                                 "be computed, use DynConnectedComponents instead.");
}

template <bool WeaklyCC>
void DynConnectedComponentsImpl<WeaklyCC>::indexEdges() {
    edgeid eid = 0;
    G->forEdges([&](node u, node v) {
        Edge edge(u, v, true);
        if (edgesMap.find(edge) == edgesMap.end()) {
            edgesMap.emplace(edge, eid);
            ++eid;
        }
    });
}

template <bool WeaklyCC>
void DynConnectedComponentsImpl<WeaklyCC>::init() {
    edgesMap.clear();
    componentPtr->reset(G->upperNodeIdBound(), none);
    isTree.assign(G->numberOfEdges(), false);
    indexEdges();
    hasRun = false;
}

template <bool WeaklyCC>
void DynConnectedComponentsImpl<WeaklyCC>::run() {
    init();

    std::queue<node> queue;
    index nComponents = 0;
    count visitedNodes = 0;
    auto &component = *componentPtr;

    for (node u : G->nodeRange()) {
        if (component[u] != none)
            continue; // u is already in some component, skip it

        // Create new component
        component.setUpperBound(nComponents + 1);
        const auto compIdx = nComponents;
        ++nComponents;

        // Start BFS and explore nodes in the new component
        queue.push(u);
        component[u] = compIdx;

        do {
            const auto curNode = queue.front();
            queue.pop();
            ++visitedNodes;

            auto visitNeighbor = [&](node, node neighbor, edgeweight, edgeid) -> void {
                if (component[neighbor] == none) {
                    queue.push(neighbor);
                    component[neighbor] = compIdx;
                    isTree[edgesMap.at(Edge(curNode, neighbor, true))] = true;
                }
            };

            G->forNeighborsOf(curNode, visitNeighbor);
            // Visit in-edges in case of weakly connected components
            if (WeaklyCC)
                G->forInNeighborsOf(curNode, visitNeighbor);

        } while (!queue.empty());

        if (visitedNodes == G->numberOfNodes())
            break; // All nodes have been visited
    }

    hasRun = true;
}

template <bool WeaklyCC>
void DynConnectedComponentsImpl<WeaklyCC>::addEdge(node u, node v) {
    edgeid eid;             // Id of the edge (u, v)
    bool isNewEdge = false; // Whether edge (u, v) was seen before

    {
        Edge edge(u, v, true);
        auto it = edgesMap.find(edge);
        if (it == edgesMap.end()) { // This edge was never seen before
            eid = static_cast<edgeid>(edgesMap.size());
            edgesMap.emplace(edge, eid);
            isNewEdge = true;
        } else
            eid = it->second;
    }

    auto &component = *componentPtr;
    // Get the IDs of u's and v's components
    index minComp, maxComp;
    std::tie(minComp, maxComp) = std::minmax(component[u], component[v]);

    // If u and v are already in the same component, we don't have to do anything
    if (maxComp == minComp) {
        if (isNewEdge)
            isTree.push_back(false);
        return;
    }

    // If u and v are in different components, we merge the two components

    // All nodes in the component with higher index (maxComp) get the smaller component id
    // (minComp). The component with highest id (lastComp) is updated to maxComp, so that we can
    // shrink partition by one.
    const index lastComp = component.upperBound() - 1;
    G->parallelForNodes([&](node w) {
        if (component[w] == maxComp)
            component[w] = minComp;
        else if (component[w] == lastComp)
            component[w] = maxComp;
    });

    // lastComp is not used anymore
    component.setUpperBound(lastComp);

    if (isNewEdge)
        isTree.push_back(true);
    else
        isTree[eid] = true;
}

template <bool WeaklyCC>
void DynConnectedComponentsImpl<WeaklyCC>::removeEdge(node u, node v) {
    const auto eid = edgesMap.at(Edge(u, v, true));

    // Edge (u, v) is not part of the spanning tree, nothing to do.
    if (!isTree[eid])
        return;

    // Edge (u, v) is part of the spanning tree, we check if an alternative path between u and v
    // exist. Mark current edge as not part of the tree.
    isTree[eid] = false;

    // Copy of the current partition. We modify a copy assuming that the component is not connected
    // anymore. We do not modify the original one because we don't know whether the current
    // partition is still connected or not.
    Partition newCmp(*componentPtr);

    // Add new component.
    const index newCmpIdx = newCmp.upperBound();

    // Ensure that the new index is not already in use.
    assert(std::all_of(componentPtr->getVector().begin(), componentPtr->getVector().end(),
                       [&](index compIdx) { return compIdx != newCmpIdx; }));

    newCmp.setUpperBound(newCmpIdx + 1);

    // Search for an alternative path from u to v
    std::queue<node> queue;
    queue.push(u);
    tmpDistances.assign(G->upperNodeIdBound(), none);
    tmpDistances[u] = 0;

    bool isConnected = false;

    do {
        const node curNode = queue.front();
        queue.pop();

        const auto curDist = tmpDistances[curNode];

        // Assign neighbor to new component
        newCmp[curNode] = newCmpIdx;

        auto visitNeighbor = [&](node neighbor) -> bool {
            if (tmpDistances[neighbor] == none) { // Unexplored neighbor
                tmpDistances[neighbor] = curDist + 1;
                if (neighbor == v) // Alternative path from u to v found
                    return true;
                queue.push(neighbor);
            }
            return false;
        };

        for (node neighbor : G->neighborRange(curNode))
            if (visitNeighbor(neighbor)) {
                isConnected = true;
                break;
            }

        if (isConnected)
            break;

        if (WeaklyCC)
            for (node neighbor : G->inNeighborRange(curNode))
                if (visitNeighbor(neighbor)) {
                    isConnected = true;
                    break;
                }
    } while (!queue.empty() && !isConnected);

    if (isConnected) // Alternative path found, update spanning tree
        updateTree(u, v);
    else // Removal disconnected the component, replace with new component
        std::swap(*componentPtr, newCmp);
}

template <bool WeaklyCC>
void DynConnectedComponentsImpl<WeaklyCC>::updateTree(node u, node v) {
    std::queue<node> q1, q2;
    q1.push(v);

    const auto vDist = tmpDistances[v];
    count level = 1;

    do {

        bool nextEdgeFound = false;
        do {
            const auto curNode = q1.front();
            q1.pop();

            auto visitNeighbor = [&](node neighbor) -> bool {
                if (tmpDistances[neighbor] == none || (vDist != tmpDistances[neighbor] + level))
                    return false; // This edge is not in the shortest path from u to v

                // Add the current edge to the spanning tree
                isTree[edgesMap.at(Edge(curNode, neighbor, true))] = true;
                nextEdgeFound = true;

                if (neighbor != u) // Target node not yet found, proceed BFS
                    q2.push(neighbor);
                return true;
            };

            for (node neighbor : G->neighborRange(curNode))
                if (visitNeighbor(neighbor))
                    break;

            if (WeaklyCC) {
                if (nextEdgeFound)
                    break;
                for (node neighbor : G->inNeighborRange(curNode))
                    if (visitNeighbor(neighbor))
                        break;
            }
        } while (!q1.empty() && !nextEdgeFound);

        // Make sure q1 is empty before swap
        std::queue<node> emptyQueue;
        std::swap(q1, emptyQueue);
        std::swap(q1, q2);

        ++level;
    } while (!q1.empty());
}

template <bool WeaklyCC>
void DynConnectedComponentsImpl<WeaklyCC>::update(GraphEvent event) {
    assureFinished();
    if (event.type == GraphEvent::EDGE_ADDITION)
        addEdge(event.u, event.v);
    else if (event.type == GraphEvent::EDGE_REMOVAL)
        removeEdge(event.u, event.v);
    else
        throw std::runtime_error("This graph event type is not supported");
}

template <bool WeaklyCC>
void DynConnectedComponentsImpl<WeaklyCC>::updateBatch(const std::vector<GraphEvent> &batch) {
    assureFinished();
    for (auto event : batch)
        update(event);
}

} // namespace DynConnectedComponentsDetails
} // namespace NetworKit

#endif /* ifndef NETWORKIT_DYN_COMPONENTS_CONNECTED_COMPONENTS_IMPL_HPP_ */
