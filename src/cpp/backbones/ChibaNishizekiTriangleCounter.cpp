/*
 * ChibaNishizekiTriangleCounter.cpp
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner
 */

#include "ChibaNishizekiTriangleCounter.h"
#include <iostream>

namespace NetworKit {

edgeCountMap ChibaNishizekiTriangleCounter::triangleCounts(const Graph& graph) {
	Graph g = graph;

	//Node attribute: marker
	std::vector<bool> nodeMarker(false);
	nodeMarker.resize(graph.numberOfNodes());

	//Edge attribute: triangle count
	edgeCountMap triangleCount;

	g.forNodes([&](node u) {
		//Mark all neighbors
		g.forNeighborsOf(u, [&](node v) {
			nodeMarker[v] = true;
		});

		//For all neighbors: check for already marked neighbors.
		g.forNeighborsOf(u, [&](node v) {
			g.forNeighborsOf(v, [&](node w) {
				if (nodeMarker[w]) {
					triangleFound(triangleCount, u, v, w);
				}
			});

			nodeMarker[v] = false;
		});

		removeNode(g, u);
	});

	return triangleCount;
}

void ChibaNishizekiTriangleCounter::removeNode(Graph& graph, const node& u) {
	//isolate the node before removing it.
	graph.forNeighborsOf(u, [&](node v) {
		graph.removeEdge(u,v);
	});

	graph.removeNode(u);
}

void ChibaNishizekiTriangleCounter::triangleFound(
		edgeCountMap& triangleCount, const node& u, const node& v,
		const node& w) {
	triangleCount[uEdge(u,v)]++;
	triangleCount[uEdge(u,w)]++;
	triangleCount[uEdge(v,w)]++;
}

/**
 * FOR TESTING/DEBUG OUTPUT PURPOSES ONLY
 */
edgeCountSet ChibaNishizekiTriangleCounter::triangleCountsDebug(const Graph& graph) {
	edgeCountMap resultMap = triangleCounts(graph);

	edgeCountSet resultSet;
	resultSet.resize(graph.numberOfEdges());

	for (auto& kv : resultMap) {
		resultSet.push_back(
				std::pair<std::pair<node, node>, count>(std::pair<node,node>(kv.first.lowNode, kv.first.highNode), kv.second));
	}

	return resultSet;
}

} /* namespace NetworKit */
