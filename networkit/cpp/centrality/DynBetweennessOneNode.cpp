/*
 * DynBetweennessOneNode.cpp
 *
 *  Created on: 10.03.2016
 *      Author: Elisabetta Bergamini
 */

#include <algorithm>
#include <queue>
#include <unordered_set>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/PrioQueue.hpp>
#include <networkit/centrality/DynBetweennessOneNode.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/SSSP.hpp>

namespace NetworKit {

DynBetweennessOneNode::DynBetweennessOneNode(Graph &G, node x) : G(G), x(x) {}

/**
 * Run method that stores a One shortest path for each node pair and stores shortest distances
 */
void DynBetweennessOneNode::run() {
    bcx = 0;
    distances.resize(G.upperNodeIdBound());
    sigma.resize(G.upperNodeIdBound());
    sigmax.resize(G.upperNodeIdBound());
    Pred.resize(G.upperNodeIdBound());

    auto computeDependencies = [&](node s) {
        // run SSSP algorithm and keep track of everything

        distances[s].resize(G.upperNodeIdBound(), infDist);
        sigma[s].resize(G.upperNodeIdBound(), 0);
        sigmax[s].resize(G.upperNodeIdBound(), 0);

        std::unique_ptr<SSSP> sssp;
        if (G.isWeighted()) {
            sssp = std::make_unique<Dijkstra>(G, s, true, true);
        } else {
            sssp = std::make_unique<BFS>(G, s, true, true);
        }

        sssp->run();

        // initialize sigma and distances for source s
        G.forNodes([&](node t) {
            distances[s][t] = sssp->distance(t);
            sssp->numberOfPaths(t).ToDouble(sigma[s][t]);
        });

        // initialize sigmax (paths that pass through x) for source s
        G.forNodes([&](node t) {
            auto paths = sssp->getPaths(t);
            for (auto path : paths) {
                if (std::find(path.begin(), path.end(), x) != path.end()) {
                    ++sigmax[s][t];
                }
            }
        });
    };

    G.forNodes(computeDependencies);

    // manually setting the sigmax to 1 for x
    sigmax[x][x] = 1;

    distancesOld = distances;

    // computing bc(x)
    bcx = 0;
    G.forNodePairs([&](node s, node y) {
        if (s != x && y != x && sigma[s][y]) {
            bcx += sigmax[s][y] / sigma[s][y];
        }
    });
    hasRun = true;
}

void DynBetweennessOneNode::update(GraphEvent event) {
    node u = event.u;
    node v = event.v;
    edgeweight weightuv = G.weight(u, v);
    if (!(event.type == GraphEvent::EDGE_ADDITION
          || (event.type == GraphEvent::EDGE_WEIGHT_INCREMENT && event.w < 0))) {
        throw std::runtime_error(
            "event type not allowed. Edge insertions and edge weight decreases only.");
    }
    if (weightuv <= distances[u][v]) {
        // initializations
        count z = G.upperNodeIdBound();
        std::vector<std::vector<node>> source_nodes(z);
        std::queue<node> Q;
        std::vector<bool> enqueued(G.upperNodeIdBound(), false);
        // queue with all visited nodes
        // if u has a new shortest path going through v, it updates the distance of u
        // and inserts u in the priority queue (or updates its priority, if already in Q)
        auto updateQueue = [&](node u) {
            if (!enqueued[u]) {
                Q.push(u);
                enqueued[u] = true;
            }
        };
        // returns smallest element in Q
        auto getMin = [&]() {
            node s = Q.front();
            Q.pop();
            return s;
        };
        // phase 1: find affected source nodes using bfs
        std::queue<node> bfsQ;
        std::vector<bool> visited(z, false);
        INFO("Phase 1. distances[", u, "][", v, "] = ", distances[u][v], ", and G.weight", u, ", ",
             v, " = ", G.weight(u, v));
        bfsQ.push(u);
        source_nodes[u].push_back(u);
        visited[u] = true;
        INFO("Entering bfs");
        while (!bfsQ.empty()) {
            node s = bfsQ.front();
            bfsQ.pop();
            DEBUG("Dequeueing node ", s);
            G.forInNeighborsOf(s, [&](node w, edgeweight) { // identify and process neighbors w of s
                if (visited[w] == false && distances[w][v] >= distances[w][u] + weightuv) {
                    bfsQ.push(w);
                    DEBUG("Pushing neighbor ", w);
                    visited[w] = true;
                    source_nodes[u].push_back(w);
                }
            });
        }
        // phase 2: for all source nodes, update distances to affected sinks
        std::vector<node> Pred(G.upperNodeIdBound());
        Pred[v] = u;
        updateQueue(v);
        INFO("Entering phase 2. source_nodes[", u, "] = ", source_nodes[u]);
        while (!Q.empty()) {
            node y = getMin();
            enqueued[y] = false;
            // update for all source nodes
            for (node s : source_nodes[Pred[y]]) {
                if (distances[s][y] >= distances[s][u] + weightuv + distances[v][y]) {
                    if (s == u && y == v) {
                        if (distances[s][y] > distances[s][u] + weightuv + distances[v][y]) {
                            // edge {u,v} is the single shortest path
                            sigma[u][v] = 1;
                            if (s == x || y == x) {
                                sigmax[u][v] = 1;
                            } else {
                                sigmax[u][v] = 0;
                            }
                        } else {
                            // edge {u,v} is new shortest path of same length
                            sigma[u][v] += 1;
                            if (s == x || y == x) {
                                sigmax[u][v] += 1;
                            } else {
                                sigmax[u][v] = 0;
                            }
                        }
                        distances[s][y] = weightuv;

                    } else {
                        if (distances[s][y] > distances[s][u] + weightuv + distances[v][y]) {
                            if (s == 1 && y == 2)
                                INFO(" GREATER s = ", s, ", y = ", y);
                            sigma[s][y] = sigma[s][u] * sigma[v][y];
                            sigmax[s][y] = sigmax[s][u] * sigma[v][y] + sigma[s][u] * sigmax[v][y];
                            distances[s][y] = distances[s][u] + weightuv + distances[v][y];
                        } else {
                            if (s == 1 && y == 2)
                                INFO(" EQUAL s = ", s, ", y = ", y);
                            sigma[s][y] += sigma[s][u] * sigma[v][y];
                            sigmax[s][y] += sigmax[s][u] * sigma[v][y] + sigma[s][u] * sigmax[v][y];
                        }
                    }
                    if (!G.isDirected()) {
                        distances[y][s] = distances[s][y];
                        sigma[y][s] = sigma[s][y];
                        sigmax[y][s] = sigmax[s][y];
                    }
                    source_nodes[y].push_back(s);
                }
                // loop over all neighbors
                G.forNeighborsOf(y, [&](node w, edgeweight weightyw) {
                    // I also check that y was a predecessor for w in the s.p. from v
                    if (distances[u][w] >= weightuv + distances[v][w]
                        && Aux::NumericTools::equal(distances[v][w], distances[v][y] + weightyw,
                                                    0.01)) {
                        Pred[w] = y;
                        updateQueue(w);
                    }
                });
            }
        }
        // computing bc(x)
        bcx = 0;
        G.forNodePairs([&](node s, node y) {
            if (s != x && y != x && sigma[s][y]) {
                bcx += sigmax[s][y] / sigma[s][y];
            }
        });
    }
}

void DynBetweennessOneNode::updateBatch(const std::vector<GraphEvent> &batch) {
    for (auto e : batch) {
        update(e);
    }
}

edgeweight DynBetweennessOneNode::getDistance(node u, node v) {
    return distances[u][v];
}

edgeweight DynBetweennessOneNode::getSigma(node u, node v) {
    return sigma[u][v];
}

edgeweight DynBetweennessOneNode::getSigmax(node u, node v) {
    return sigmax[u][v];
}

edgeweight DynBetweennessOneNode::getbcx() {
    return bcx;
}

} /* namespace NetworKit */
