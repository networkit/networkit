/*
 * DynDijkstra.cpp
 *
 *  Created on: 21.07.2014
 *      Author: ebergamini
 */

#include <queue>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/DynDijkstra.hpp>

namespace NetworKit {

DynDijkstra::DynDijkstra(const Graph &G, node source, bool storePredecessors)
    : DynSSSP(G, source, storePredecessors), color(G.upperNodeIdBound(), WHITE),
      heap(Aux::LessInVector<edgeweight>(distances)), updateDistances(G.upperNodeIdBound()),
      updateHeap(Aux::LessInVector<edgeweight>(updateDistances)) {}

void DynDijkstra::run() {

    color.clear();
    count n = G->upperNodeIdBound();
    color.resize(n, WHITE);

    // init distances
    distances.clear();
    std::fill(distances.begin(), distances.end(), infDist);

    if (distances.size() < G->upperNodeIdBound())
        distances.resize(G->upperNodeIdBound(), infDist);

    sumDist = 0.;
    reachedNodes = 1;

    if (storePreds) {
        previous.clear();
        previous.resize(G->upperNodeIdBound());
    }

    npaths.clear();
    npaths.resize(G->upperNodeIdBound(), 0);
    npaths[source] = 1;

    // priority queue with distance-node pairs
    distances[source] = 0.;
    heap.clear();
    heap.reserve(n);
    heap.push(source);

    auto initPath = [&](node u, node v) {
        if (storePreds) {
            previous[v] = {u};
        }
        npaths[v] = npaths[u];
    };

    TRACE("traversing graph");
    do {
        TRACE("pq size: ", heap.size());
        node u = heap.extract_top();
        sumDist += distances[u];
        TRACE("current node in Dijkstra: ", u);

        G->forNeighborsOf(u, [&](node v, edgeweight w) {
            double newDist = distances[u] + w;
            if (distances[v] == infDist) {
                distances[v] = newDist;
                heap.push(v);
                ++reachedNodes;
                initPath(u, v);
            } else if (distances[v] > newDist) {
                initPath(u, v);
                distances[v] = newDist;
                heap.update(v);
            } else if (Aux::NumericTools::logically_equal(distances[v], newDist)) {
                if (storePreds)
                    previous[v].push_back(u);
                npaths[v] += npaths[u];
            }
        });
    } while (!heap.empty());

    hasRun = true;
}

void DynDijkstra::update(GraphEvent e) {
    std::vector<GraphEvent> batch(1);
    batch[0] = e;
    updateBatch(batch);
}

void DynDijkstra::updateBatch(const std::vector<GraphEvent> &batch) {
    mod = false;
    // priority queue with distance-node pairs
    updateHeap.clear();
    updateHeap.reserve(G->upperNodeIdBound());

    // queue with all visited nodes
    std::queue<node> visited;
    // if u has a new shortest path going through v, it updates the distance of u
    // and inserts u in the priority queue (or updates its priority, if already in Q)
    auto updateQueue = [&](node v, edgeweight dist) {
        if (color[v] == WHITE) {
            updateDistances[v] = dist;
            updateHeap.push(v);
            color[v] = GRAY;
        } else {
            updateDistances[v] = dist;
            updateHeap.update(v);
        }
    };

    // initialization
    for (GraphEvent edge : batch) {

        if (!G->hasNode(edge.u) || !G->hasNode(edge.v)) {
            throw std::runtime_error("Graph update not allowed or invalid nodes");
        }
        if (distances[edge.u] == infDist && distances[edge.v] == infDist)
            continue;

        node vUp, vDown;

        if (distances[edge.u] < distances[edge.v]) {
            vUp = edge.u;
            vDown = edge.v;
        } else if (distances[edge.u] >= distances[edge.v]) {
            vUp = edge.v;
            vDown = edge.u;
        } else {
            continue;
        }

        if (edge.type == GraphEvent::EDGE_REMOVAL) {
            if (distances[vDown] == infDist)
                continue;
            updateQueue(vDown, distances[vDown]);
            continue;
        } else if (distances[vUp] + edge.w < distances[vDown]) {
            updateQueue(vDown, distances[vUp] + edge.w);
        } else {
            updateQueue(vDown, distances[vDown]);
        }
    }

    while (!updateHeap.empty()) {
        mod = true;
        node current = updateHeap.extract_top();
        // the previous vector for the current node is potentially outdated and will be refilled
        // during the update
        if (storePreds) {
            previous[current].clear();
        }
        if (color[current] == BLACK
            && Aux::NumericTools::logically_equal(updateDistances[current], distances[current])) {
            auto tmp = npaths[current];
            npaths[current] = 0;
            G->forInNeighborsOf(current, [&](node current, node z, edgeweight w) {
                // if z is a predecessor for current update the shortest paths
                if (Aux::NumericTools::logically_equal(distances[current], distances[z] + w)) {
                    if (storePreds) {
                        previous[current].push_back(z);
                    }
                    npaths[current] += npaths[z];
                }
            });
            if (tmp != npaths[current]) {
                G->forNeighborsOf(current, [&](node z, edgeweight w) {
                    // current was a predecessor for z
                    if (Aux::NumericTools::ge(distances[current] + w, distances[z], distEpsilon))
                        updateQueue(z, distances[current] + w);
                });
            }
            continue;
        }
        edgeweight k = updateDistances[current];
        edgeweight con = infDist;
        // check whether there are still predecessors for current node
        G->forInNeighborsOf(current, [&](node z, edgeweight w) {
            if (distances[z] != infDist && distances[z] + w < con) {
                con = distances[z] + w;
            }
        });
        npaths[current] = 0;
        if (Aux::NumericTools::logically_equal(con, k)) {

            distances[current] = k;
            visited.push(current);
            color[current] = BLACK;
            G->forInNeighborsOf(current, [&](node current, node z, edgeweight w) {
                // if z is a predecessor for current update the shortest paths
                if (Aux::NumericTools::logically_equal(distances[current], distances[z] + w)) {
                    if (storePreds) {
                        previous[current].push_back(z);
                    }
                    npaths[current] += npaths[z];
                }
            });

            G->forNeighborsOf(current, [&](node z, edgeweight w) {
                // current is a predecessor for z
                if (distances[current] == infDist
                    || Aux::NumericTools::ge(distances[z], k + w, distEpsilon)) {
                    updateQueue(z, k + w);
                }
            });

        } else { // con > k
            if (distances[current] != infDist) {
                distances[current] = infDist;
                G->forNeighborsOf(current, [&](node z, edgeweight w) {
                    // current was a predecessor for z
                    if (Aux::NumericTools::ge(distances[z], k + w, distEpsilon)) {
                        updateQueue(z, k + w);
                    }
                });
            }
            if (con != infDist) {
                updateQueue(current, con);
            }
        }
    }

    // reset colors
    while (!visited.empty()) {
        node w = visited.front();
        visited.pop();
        color[w] = WHITE;
    }
}

constexpr edgeweight DynDijkstra::infDist;
constexpr edgeweight DynDijkstra::distEpsilon;

} /* namespace NetworKit */
