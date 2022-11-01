/*
 * DynBFS.cpp
 *
 *  Created on: 17.07.2014
 *      Author: cls, ebergamini
 */

#include <queue>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/DynBFS.hpp>

namespace NetworKit {

DynBFS::DynBFS(const Graph &G, node s, bool storePredecessors)
    : DynSSSP(G, s, storePredecessors), color(G.upperNodeIdBound(), WHITE) {}

void DynBFS::run() {
    count n = G->upperNodeIdBound();
    distances.clear();
    distances.resize(n, infDist);
    color.clear();
    color.resize(n, WHITE);
    std::vector<bool> visited;
    visited.resize(n, false);

    if (storePreds) {
        previous.clear();
        previous.resize(n);
    }

    npaths.clear();
    npaths.resize(n, 0);
    npaths[source] = 1;
    std::queue<node> q;
    q.push(source);
    visited[source] = true;
    distances[source] = 0.0;
    maxDistance = 0;
    while (!q.empty()) {
        node u = q.front();
        q.pop();
        if (distances[u] > maxDistance)
            maxDistance = distances[u];
        // insert untouched neighbors into queue
        G->forNeighborsOf(u, [&](node v) {
            if (distances[v] == infDist) {
                q.push(v);
                distances[v] = distances[u] + 1.;
                sumDist += distances[v];
                ++reachedNodes;
                if (storePreds) {
                    previous[v] = {u};
                }
                npaths[v] = npaths[u];
            } else if (distances[v] == distances[u] + 1.) {
                if (storePreds)
                    previous[v].push_back(u); // additional predecessor
                npaths[v] += npaths[u];       // all the shortest paths to u are also
                                              // shortest paths to v now
            }
            if (!visited[v]) {
                visited[v] = true;
            }
        });
    }
    hasRun = true;
}

void DynBFS::update(GraphEvent e) {
    std::vector<GraphEvent> batch(1);
    batch[0] = e;
    updateBatch(batch);
}

void DynBFS::updateBatch(const std::vector<GraphEvent> &batch) {
    std::vector<std::queue<node>> queues(maxDistance + 1);
    mod = false;
    // insert nodes from the batch whose distance has changed (affected nodes) into the queues
    for (GraphEvent edge : batch) {
        if (edge.type == GraphEvent::EDGE_ADDITION) {
            if (distances[edge.u] == infDist && (G->isDirected() || distances[edge.v] == infDist))
                continue;
            if (!G->isDirected() && distances[edge.u] > distances[edge.v]) {
                queues[distances[edge.v] + 1].push(edge.u);
            }
            if (distances[edge.v] > distances[edge.u]) {
                queues[distances[edge.u] + 1].push(edge.v);
            }
        } else if (edge.type == GraphEvent::EDGE_REMOVAL) { // edge deletion.
            if (distances[edge.u] == infDist || distances[edge.v] == infDist)
                continue;
            if (!G->isDirected() && distances[edge.u] > distances[edge.v]) {
                queues[distances[edge.v] + 1].push(edge.u);
            }
            if (distances[edge.v] > distances[edge.u]) {
                queues[distances[edge.u] + 1].push(edge.v);
            }
        } else
            throw std::runtime_error(
                "Graph update not allowed: only edge insertions and edge deletions");
    }
    // extract nodes from the queues and scan incident edges
    std::queue<node> visited;

    for (count m = 1; m < maxDistance; m++) {
        if (m >= maxDistance - 1 && (!queues[m].empty() || !queues[m + 1].empty())) {
            maxDistance++;
            queues.emplace_back();
        }
        mod = mod || (!queues[m].empty());
        while (!queues[m].empty()) {
            node w = queues[m].front();
            queues[m].pop();
            if (color[w] == BLACK) {
                continue;
            }
            edgeweight con = infDist;
            // check whether there are still predecessors for w at level m
            G->forInNeighborsOf(w, [&](node z) {
                if (distances[z] != infDist && distances[z] + 1 < con) {
                    con = distances[z] + 1;
                }
            });
            if (con == m) {
                npaths[w] = 0;
                distances[w] = m;
                visited.push(w);
                color[w] = BLACK;
                G->forInNeighborsOf(w, [&](node z) {
                    // z is a predecessor for w
                    if (distances[w] == distances[z] + 1)
                        npaths[w] += npaths[z];
                });
                G->forNeighborsOf(w, [&](node z) {
                    // w is a predecessor for z
                    if (distances[w] == infDist || distances[z] >= double(m) + 1.0) {
                        queues[m + 1].push(z);
                        color[z] = GRAY;
                    }
                });
            } else {
                assert(con > m);
                if (distances[w] != infDist) {
                    npaths[w] = 0;
                    distances[w] = infDist;
                    G->forNeighborsOf(w, [&](node z) {
                        // w was a predecessor for z
                        if (distances[z] > m && distances[z] < infDist) {
                            queues[m + 1].push(z);
                            color[z] = GRAY;
                        }
                    });
                }
                if (con != infDist) {
                    if (con > maxDistance) {
                        for (count i = maxDistance; i < con; i++)
                            queues.emplace_back();
                        maxDistance = con;
                    }
                    queues[con].push(w);
                }
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

constexpr edgeweight DynBFS::infDist;

} /* namespace NetworKit */
