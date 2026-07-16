/*
 * BFSImpl.hpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_DISTANCE_BFS_IMPL_HPP_
#define NETWORKIT_DISTANCE_BFS_IMPL_HPP_

#include <queue>

namespace NetworKit {

template <class GraphT>
BreadthFirstSearch<GraphT>::BreadthFirstSearch(const GraphT &G, NodeT source, bool storePaths,
                                               bool storeNodesSortedByDistance, NodeT target)
    : SingleSourceShortestPaths<GraphT>(G, source, storePaths, storeNodesSortedByDistance, target) {
}

template <class GraphT>
void BreadthFirstSearch<GraphT>::run() {
    count z = this->G->upperNodeIdBound();
    this->reachedNodes = 1;
    this->sumDist = 0.;

    const auto infDist = std::numeric_limits<EdgeWeightT>::max();
    std::fill(this->distances.begin(), this->distances.end(), infDist);

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
    this->distances[this->source] = 0.;

    bool breakWhenFound = (this->target != nullNodeId);
    while (!q.empty()) {
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
                this->distances[v] = this->distances[u] + 1.;
                this->sumDist += this->distances[v];
                ++this->reachedNodes;
                if (this->storePaths) {
                    this->previous[v] = {u};
                    this->npaths[v] = this->npaths[u];
                }
            } else if (this->storePaths && (this->distances[v] == this->distances[u] + 1.)) {
                // additional predecessor
                this->previous[v].push_back(u);
                // all the shortest paths to u are also shortest paths to v now
                this->npaths[v] += this->npaths[u];
            }
        });
    }

    this->hasRun = true;
}

} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_BFS_IMPL_HPP_
