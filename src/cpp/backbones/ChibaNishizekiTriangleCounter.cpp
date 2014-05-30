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

/**
 * FOR TESTING PURPOSES ONLY
 */
edgeCountSet ChibaNishizekiTriangleCounter::triangleCountsDebug(const Graph& graph) {
	edgeCountMap resultMap = triangleCounts(graph);

	edgeCountSet resultSet;
	resultSet.resize(graph.numberOfEdges());

	for (auto& kv : resultMap) {
		resultSet.push_back(
				std::pair<std::pair<node, node>, count>(kv.first, kv.second));
	}

	return resultSet;
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
	//TODO: remove this ugliness. maybe custom "edge" key class?
	if (u < v)
		triangleCount[std::pair<node,node>(u,v)]++;
	else
		triangleCount[std::pair<node,node>(v,u)]++;

	if (u < w)
		triangleCount[std::pair<node,node>(u,w)]++;
	else
		triangleCount[std::pair<node,node>(w,u)]++;

	if (v < w)
		triangleCount[std::pair<node,node>(v,w)]++;
	else
		triangleCount[std::pair<node,node>(w,v)]++;
}

} /* namespace NetworKit */
