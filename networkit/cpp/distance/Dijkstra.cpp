/*
 * Dijkstra.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning, Christian Staudt
 */

#include <algorithm>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/distance/Dijkstra.hpp>

namespace NetworKit {

Dijkstra::Dijkstra(const Graph &G, node source, bool storePaths,
                   bool storeNodesSortedByDistance, node target)
    : SSSP(G, source, storePaths, storeNodesSortedByDistance, target),
      heap(Compare(distances)) {}

void Dijkstra::run() {

    TRACE("initializing Dijkstra data structures");
    // init distances
    if (distances.size() < G->upperNodeIdBound()) {
        distances.resize(G->upperNodeIdBound(),
                         std::numeric_limits<double>::max());
        visited.resize(G->upperNodeIdBound(), ts);
    }

    sumDist = 0.;
    reachedNodes = 1;

    if (ts++ == 255) {
        ts = 1;
        std::fill(visited.begin(), visited.end(), 0);
    }

    if (storePaths) {
        previous.clear();
        previous.resize(G->upperNodeIdBound());
        npaths.clear();
        npaths.resize(G->upperNodeIdBound(), 0);
        npaths[source] = 1;
    }

    if (storeNodesSortedByDistance) {
        nodesSortedByDistance.clear();
        nodesSortedByDistance.reserve(G->upperNodeIdBound());
    }

    // priority queue with distance-node pairs
    distances[source] = 0.;
    visited[source] = ts;
    heap.push(source);

    auto initPath = [&](node u, node v) {
        if (storePaths) {
            previous[v] = {u};
            npaths[v] = npaths[u];
        }
    };
    bool breakWhenFound = (target != none);
    TRACE("traversing graph");
    do {
        TRACE("pq size: ", heap.size());
        node u = heap.extract_top();
        sumDist += distances[u];
        TRACE("current node in Dijkstra: ", u);
        TRACE("pq size: ", heap.size());
        if ((breakWhenFound && target == u) || visited[u] != ts)
            break;

        if (storeNodesSortedByDistance)
            nodesSortedByDistance.push_back(u);

        G->forNeighborsOf(u, [&](node v, edgeweight w) {
            double newDist = distances[u] + w;
            if (ts != visited[v]) {
                visited[v] = ts;
                distances[v] = newDist;
                heap.push(v);
                ++reachedNodes;
                if (storePaths)
                    initPath(u, v);
            } else if (distances[v] > newDist) {
                if (storePaths)
                    initPath(u, v);
                distances[v] = newDist;
                heap.update(v);
            } else if (storePaths && distances[v] == newDist) {
                previous[v].push_back(u);
                npaths[v] += npaths[u];
            }
        });
    } while (!heap.empty());

    hasRun = true;
}
} // namespace NetworKit
