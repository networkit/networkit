/*
 * BFS.hpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_DISTANCE_BFS_HPP_
#define NETWORKIT_DISTANCE_BFS_HPP_

#include <algorithm>
#include <queue>

#include <networkit/distance/SSSP.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * The BFS class is used to do a breadth-first search on a Graph from a given
 * source node.
 */
template <class GraphT>
class BreadthFirstSearch final : public SingleSourceShortestPaths<GraphT> {
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;
    static constexpr NodeT nullNodeId = NullNodeId<NodeT>;

public:
    /**
     * Constructs the BreadthFirstSearch class for @a G and source node @a source.
     *
     * @param G The graph
     * @param source The source node of the breadth-first search
     * @param storePaths Paths are reconstructable and the number of paths is
     * stored.
     * @param storeNodesSortedByDistance Store a vector of nodes ordered in
     * increasing distance from the source.
     * @param target The target node.
     */
    BreadthFirstSearch(const GraphT &G, NodeT source, bool storePaths = true,
                       bool storeNodesSortedByDistance = false, NodeT target = nullNodeId)
        : SingleSourceShortestPaths<GraphT>(G, source, storePaths, storeNodesSortedByDistance,
                                            target) {}

    /**
     * Breadth-first search from @a source.
     */
    void run() override;
};

template <class GraphT>
void BreadthFirstSearch<GraphT>::run() {
    constexpr EdgeWeightT weightOne{1};
    const count z = this->G->upperNodeIdBound();
    this->reachedNodes = 1;
    this->sumDist = 0.;

    const auto infDist = std::numeric_limits<EdgeWeightT>::max();
    std::ranges::fill(this->distances, infDist);

    if (this->distances.size() < z)
        this->distances.resize(z, infDist);

    if (this->storePaths) {
        this->previous.clear();
        this->previous.resize(z);
        this->npaths.clear();
        this->npaths.resize(z, 0);
        this->npaths[this->source] = 1;
    }

    if (this->storeNodesSortedByDistance) {
        std::vector<NodeT> empty;
        std::swap(this->nodesSortedByDistance, empty);
    }

    std::queue<NodeT> q;
    q.push(this->source);
    this->distances[this->source] = EdgeWeightT{0};

    bool breakWhenFound = (this->target != nullNodeId);
    do {
        NodeT u = q.front();
        q.pop();

        if (this->storeNodesSortedByDistance) {
            this->nodesSortedByDistance.push_back(u);
        }
        if (breakWhenFound && this->target == u) {
            break;
        }

        // insert untouched neighbors into queue
        this->G->forNeighborsOf(u, [&](NodeT v) {
            if (this->distances[v] == infDist) {
                q.push(v);
                this->distances[v] = this->distances[u] + weightOne;
                this->sumDist += this->distances[v];
                ++this->reachedNodes;
                if (this->storePaths) {
                    this->previous[v] = {u};
                    this->npaths[v] = this->npaths[u];
                }
            } else if (this->storePaths && (this->distances[v] == this->distances[u] + weightOne)) {
                // additional predecessor
                this->previous[v].push_back(u);
                // all the shortest paths to u are also shortest paths to v now
                this->npaths[v] += this->npaths[u];
            }
        });
    } while (!q.empty());

    this->hasRun = true;
}

using BFS = BreadthFirstSearch<AdjListGraph<node, edgeweight>>;

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_BFS_HPP_
