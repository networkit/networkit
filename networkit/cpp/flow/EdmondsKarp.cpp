/*
 * EdmondsKarp.cpp
 *
 *  Created on: 11.06.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "EdmondsKarp.h"

namespace NetworKit {

edgeweight EdmondsKarp::BFS(const Graph &graph, std::vector<edgeweight> &flow, std::vector<edgeweight> &residFlow, node source, node sink, std::vector<node> &pred) const {
	pred = std::vector<node>(graph.upperNodeIdBound(), none);
	std::vector<edgeweight> gain(graph.upperNodeIdBound(), 0);

	std::queue<node> Q;
	Q.push(source);
	pred[source] = source;
	gain[source] = std::numeric_limits<edgeweight>::max();
	while (!Q.empty()) {
		node u = Q.front(); Q.pop();

		bool sinkReached = false;
		graph.forNeighborsOf(u, [&](node _u, node v, edgeweight weight, edgeid eid){
			if ((
			(u >= v && flow[eid] < weight) || (u < v && residFlow[eid] < weight)
			)&& pred[v] == none) { // only add those neighbors with rest capacity and which were not discovered yet
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
	}

	return 0.0;
}

edgeweight EdmondsKarp::solveMaxFlow(const Graph &graph, const node source, const node sink, std::vector<edgeweight> &flow) const {
	if (!graph.hasEdgeIds()) { throw std::runtime_error("edges have not been indexed - call indexEdges first"); }
	flow.clear();
	flow.resize(graph.upperEdgeIdBound(), 0.0);

	std::vector<edgeweight> residFlow(graph.upperEdgeIdBound(), 0.0);

	edgeweight maxFlow = 0;
	while (true) {
		std::vector<node> pred;
		edgeweight gain = BFS(graph, flow, residFlow, source, sink, pred);
		if (gain == 0) break;

		maxFlow += gain;
		node v = sink;
		while (v != source) {
			node u = pred[v];
			edgeid eid = graph.edgeId(u, v);
			if (u >= v) {
				flow[eid] += gain;
				residFlow[eid] -= gain;
			} else {
				flow[eid] -= gain;
				residFlow[eid] += gain;
			}
			v = u;
		}
	}

	graph.parallelForEdges([&](node u, node v, edgeid eid) {
		flow[eid] = std::max(flow[eid], residFlow[eid]);
	});

	return maxFlow;
}

void EdmondsKarp::computeSourceSet(const Graph &graph, const node source, const node sink, const std::vector<edgeweight> &flow, std::vector<node> &sourceSet) const {
	// perform bfs from source
	std::vector<bool> visited(graph.upperNodeIdBound(), false);
	sourceSet.clear();
	std::queue<node> Q;
	Q.push(source);
	visited[source] = true;
	while (!Q.empty()) {
		node u = Q.front(); Q.pop();
		sourceSet.push_back(u);

		graph.forNeighborsOf(u, [&](node _u, node v, edgeweight weight, edgeid eid) {
			if (!visited[v] && flow[eid] < weight) {
				Q.push(v);
				visited[v] = true;
			}
		});
	}
}

edgeweight EdmondsKarp::run(const Graph &graph, const node source, const node sink) const {
	std::vector<edgeweight> flow;
	return solveMaxFlow(graph, source, sink, flow);
}

edgeweight EdmondsKarp::run(const Graph &graph, const node source, const node sink, std::vector<node> &sourceSet) const {
	std::vector<edgeweight> flow;
	edgeweight maxFlow = solveMaxFlow(graph, source, sink, flow);
	computeSourceSet(graph, source, sink, flow, sourceSet);

	return maxFlow;
}

edgeweight EdmondsKarp::run(const NetworKit::Graph &graph, const node source, const node sink, std::vector< NetworKit::edgeweight > &flow) const {
	edgeweight maxFlow = solveMaxFlow(graph, source, sink, flow);

	return maxFlow;
}

edgeweight EdmondsKarp::run(const NetworKit::Graph &graph, node source, node sink, std::vector< NetworKit::node > &sourceSet, std::vector< NetworKit::edgeweight > &flow) const {
	edgeweight maxFlow = solveMaxFlow(graph, source, sink, flow);

	computeSourceSet(graph, source, sink, flow, sourceSet);

	return maxFlow;
}


} /* namespace NetworKit */
