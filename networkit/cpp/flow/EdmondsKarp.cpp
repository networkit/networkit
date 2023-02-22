/*
 * EdmondsKarp.cpp
 *
 *  Created on: 11.06.2014
 *     Authors: Michael Wegner <michael.wegner@student.kit.edu>
 *              Michael Hamann <michael.hamann@kit.edu>
 */

#include <algorithm>
#include <limits>
#include <queue>
#include <stdexcept>

#include <networkit/flow/EdmondsKarp.hpp>

namespace NetworKit {

EdmondsKarp::EdmondsKarp(const Graph &graph, node source, node sink)
    : graph(&graph), source(source), sink(sink) {}

edgeweight EdmondsKarp::BFS(const std::vector<edgeweight> &residFlow,
                            std::vector<node> &pred) const {
    std::fill(pred.begin(), pred.end(), none);
    pred.resize(graph->upperNodeIdBound(), none);
    std::vector<edgeweight> gain(graph->upperNodeIdBound(), 0);

    std::queue<node> Q;
    Q.push(source);
    pred[source] = source;
    gain[source] = std::numeric_limits<edgeweight>::max();
    do {
        node u = Q.front();
        Q.pop();

        bool sinkReached = false;
        graph->forNeighborsOf(u, [&](node, node v, edgeweight weight, edgeid eid) {
            if (((u >= v && flow[eid] < weight) || (u < v && residFlow[eid] < weight))
                && pred[v] == none) {
                // only add those neighbors with rest capacity and which were not discovered yet
                pred[v] = u;
                gain[v] = std::min(gain[u], weight - (u >= v ? flow[eid] : residFlow[eid]));

                if (v != sink && !sinkReached) {
                    Q.push(v);
                } else {
                    sinkReached = true;
                }
            }
        });

        if (sinkReached) {
            return gain[sink];
        }
    } while (!Q.empty());

    return 0.0;
}

edgeweight EdmondsKarp::directedBFS(const std::vector<count> &reverseEdges,
                                    std::vector<node> &pred) const {
    std::fill(pred.begin(), pred.end(), none);
    pred.resize(graph->upperNodeIdBound(), none);
    std::vector<edgeweight> gain(graph->upperNodeIdBound(), std::numeric_limits<edgeweight>::max());

    std::queue<node> Q;
    Q.push(source);
    pred[source] = source;
    gain[source] = std::numeric_limits<edgeweight>::max();;

    while (!Q.empty()) {
        node u = Q.front();
        Q.pop();
        bool sinkReached = false;

        graph->forNeighborsOf(u, [&](node, node v, edgeweight weight, edgeid eid) {
            if (!sinkReached && pred[v] == none && flow[eid] < weight) {
                pred[v] = u;
                if (v != sink && !sinkReached) {
                    Q.push(v);
                } else {
                    sinkReached = true;
                }

                auto residCapa =
                    weight - flow[eid] + (reverseEdges[eid] == none ? 0 : flow[reverseEdges[eid]]);
                gain[v] = std::min(gain[u], residCapa);
            }
        });

        if (sinkReached) {
            return gain[sink];
        }

        graph->forInNeighborsOf(u, [&](node, node v, edgeweight weight, edgeid eid) {
            if (!sinkReached && (reverseEdges[eid] == none) && (flow[eid] > 0)
                && (pred[v] == none)) {
                pred[v] = u;
                gain[v] = std::min(gain[u], flow[eid]);

                if (v != sink && !sinkReached) {
                    Q.push(v);
                } else {
                    sinkReached = true;
                }
            }
        });

        if (sinkReached) {
            return gain[sink];
        }
    }

    return 0.0;
}

void EdmondsKarp::run() {
    if (!graph->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }
    flow.clear();
    flow.resize(graph->upperEdgeIdBound(), 0.0);
    flowValue = 0;

    if (graph->isDirected()) {
        runDirected();
    } else {
        runUndirected();
    }

    hasRun = true;
}

void EdmondsKarp::runDirected() {
    std::vector<count> reverseEdges(graph->upperEdgeIdBound(), none);
    graph->parallelForEdges([&](node u, node v, edgeid eid) {
        reverseEdges[eid] = graph->hasEdge(v, u) ? graph->edgeId(v, u) : none;
    });
    std::vector<node> pred;
    edgeweight gain = directedBFS(reverseEdges, pred);

    while (gain > 0) {
        flowValue += gain;
        node v = sink;
        while (v != source) {
            node u = pred[v];

            if (graph->hasEdge(u, v)) {
                edgeid eid = graph->edgeId(u, v);
                auto reverseId = reverseEdges[eid];
                if (reverseId != none) {
                    auto reverseFlow = flow[reverseId];
                    assert(!((reverseFlow > 0) && (flow[eid] > 0)));
                    flow[reverseId] = std::max(reverseFlow - gain, 0.);
                    flow[eid] += gain - reverseFlow;
                } else {
                    flow[eid] += gain;
                    assert(flow[eid] <= graph->weight(u, v));
                }
            } else {
                edgeid eid = graph->edgeId(v, u);
                flow[eid] -= gain;
                assert(flow[eid] >= 0);
            }

            v = u;
        }
        gain = directedBFS(reverseEdges, pred);
    }
}

void EdmondsKarp::runUndirected() {
    std::vector<edgeweight> residFlow(graph->upperEdgeIdBound(), 0.0);

    do {
        std::vector<node> pred;
        edgeweight gain = BFS(residFlow, pred);
        if (gain == 0)
            break;

        flowValue += gain;
        node v = sink;
        while (v != source) {
            node u = pred[v];
            edgeid eid = graph->edgeId(u, v);
            if (u >= v) {
                flow[eid] += gain;
                residFlow[eid] -= gain;
            } else {
                flow[eid] -= gain;
                residFlow[eid] += gain;
            }
            v = u;
        }
    } while (true);

    graph->parallelForEdges(
        [&](node, node, edgeid eid) { flow[eid] = std::max(flow[eid], residFlow[eid]); });
}

edgeweight EdmondsKarp::getMaxFlow() const {
    assureFinished();
    return flowValue;
}

std::vector<node> EdmondsKarp::getSourceSet() const {
    assureFinished();
    // perform bfs from source
    std::vector<bool> visited(graph->upperNodeIdBound(), false);
    std::vector<node> sourceSet;

    std::queue<node> Q;
    Q.push(source);
    visited[source] = true;
    while (!Q.empty()) {
        node u = Q.front();
        Q.pop();
        sourceSet.push_back(u);

        graph->forNeighborsOf(u, [&](node, node v, edgeweight weight, edgeid eid) {
            if (!visited[v] && flow[eid] < weight) {
                Q.push(v);
                visited[v] = true;
            }
        });
        if (graph->isDirected()) {
            // Follow backwards edges in residual network
            graph->forInNeighborsOf(u, [&](node, node v, edgeid eid) {
                if (!visited[v] && flow[eid] > 0) {
                    Q.push(v);
                    visited[v] = true;
                }
            });
        }
    }

    return sourceSet;
}

edgeweight EdmondsKarp::getFlow(node u, node v) const {
    assureFinished();
    return flow[graph->edgeId(u, v)];
}

const std::vector<edgeweight> &EdmondsKarp::getFlowVector() const {
    assureFinished();
    return flow;
}

} /* namespace NetworKit */
