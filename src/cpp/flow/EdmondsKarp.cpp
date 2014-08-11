/*
 * EdmondsKarp.cpp
 *
 *  Created on: 11.06.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "EdmondsKarp.h"

namespace NetworKit {

index EdmondsKarp::getEdgeIdx(const Graph &graph, node u, node v) const {
	index idx = 0;
	bool found = false;
	graph.forNeighborsOf(u, [&](node neighbor){
		if (neighbor == v) {
			found = true;
		}

		if (!found) {
			idx++;
		}
	});

	return idx;
}

edgeweight EdmondsKarp::BFS(const Graph &graph, std::vector<std::vector<edgeweight>> &flow, node source, node sink, std::vector<node> &pred) const {
	pred = std::vector<node>(graph.numberOfNodes(), none);
	std::vector<edgeweight> gain(graph.numberOfNodes(), 0);

	std::queue<node> Q;
	Q.push(source);
	pred[source] = source;
	gain[source] = std::numeric_limits<edgeweight>::max();
	while (!Q.empty()) {
		node u = Q.front(); Q.pop();

		bool sinkReached = false;
		graph.forNeighborsOf(u, [&](node v){
			index edgeIdx = getEdgeIdx(graph, u, v);
			if (flow[u][edgeIdx] < graph.weight(u,v) && pred[v] == none) { // only add those neighbors with rest capacity and which were not discovered yet
				pred[v] = u;
				gain[v] = std::min(gain[u], graph.weight(u,v) - flow[u][edgeIdx]);

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

edgeweight EdmondsKarp::solveMaxFlow(const Graph &graph, const node source, const node sink, std::vector<std::vector<edgeweight>> &flow) const {
	flow = std::vector<std::vector<edgeweight>>(graph.numberOfNodes());
	graph.forNodes([&](node u){
		flow[u] = std::vector<edgeweight>(graph.degree(u), 0.0);
	});

	edgeweight maxFlow = 0;
	while (true) {
		std::vector<node> pred;
		edgeweight gain = BFS(graph, flow, source, sink, pred);
		if (gain == 0) break;

		maxFlow += gain;
		node v = sink;
		while (v != source) {
			node u = pred[v];
			flow[u][getEdgeIdx(graph, u, v)] += gain;
			flow[v][getEdgeIdx(graph, v, u)] += gain; // undirected!
			v = u;
		}
	}

	return maxFlow;
}

void EdmondsKarp::computeSourceSet(const Graph &graph, const node source, const node sink, const std::vector<std::vector<edgeweight>> &flow, std::vector<node> &sourceSet) const {
	// perform bfs from source
	std::vector<bool> visited(graph.numberOfNodes(), false);
	sourceSet.clear();
	std::queue<node> Q;
	Q.push(source);
	visited[source] = true;
	while (!Q.empty()) {
		node u = Q.front(); Q.pop();
		sourceSet.push_back(u);

		graph.forNeighborsOf(u, [&](node v) {
			if (!visited[v] && flow[u][getEdgeIdx(graph, u, v)] < graph.weight(u,v)) {
				Q.push(v);
				visited[v] = true;
			}
		});
	}
}

edgeweight EdmondsKarp::run(const Graph &graph, const node source, const node sink) const {
	std::vector<std::vector<edgeweight>> flow;
	return solveMaxFlow(graph, source, sink, flow);
}

edgeweight EdmondsKarp::run(const Graph &graph, const node source, const node sink, std::vector<node> &sourceSet) const {
	std::vector<std::vector<edgeweight>> flow;
	edgeweight maxFlow = solveMaxFlow(graph, source, sink, flow);
	computeSourceSet(graph, source, sink, flow, sourceSet);

	return maxFlow;
}

edgeweight EdmondsKarp::run(Graph &graph, const node source, const node sink, int &attribute_id) const {
	std::vector<std::vector<edgeweight>> flow;
	edgeweight maxFlow = solveMaxFlow(graph, source, sink, flow);

	attribute_id = graph.addEdgeAttribute_double(0.0);
	graph.forEdges([&](node u, node v) {
		graph.setAttribute_double(u, v, attribute_id, flow[u][getEdgeIdx(graph, u, v)]);
	});

	return maxFlow;
}

edgeweight EdmondsKarp::run(Graph &graph, node source, node sink, std::vector<node> &sourceSet, int &attribute_id) const {
	std::vector<std::vector<edgeweight>> flow;
	edgeweight maxFlow = solveMaxFlow(graph, source, sink, flow);

	attribute_id = graph.addEdgeAttribute_double(0.0);
	graph.forEdges([&](node u, node v) {
		graph.setAttribute_double(u, v, attribute_id, flow[u][getEdgeIdx(graph, u, v)]);
	});

	computeSourceSet(graph, source, sink, flow, sourceSet);

	return maxFlow;
}


} /* namespace NetworKit */
