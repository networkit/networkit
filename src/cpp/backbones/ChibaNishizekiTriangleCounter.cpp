/*
 * ChibaNishizekiTriangleCounter.cpp
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner
 */

#include "ChibaNishizekiTriangleCounter.h"

namespace NetworKit {

ChibaNishizekiTriangleCounter::ChibaNishizekiTriangleCounter() {
}

ChibaNishizekiTriangleCounter::~ChibaNishizekiTriangleCounter() {
}

edgeCountSet ChibaNishizekiTriangleCounter::triangleCounts(const Graph& graph) {
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
				triangleCount[std::pair<node,node>(u,v)]++;
				//triangleCount[std::pair<node,node>(v,u)]++;
				triangleCount[std::pair<node,node>(u,w)]++;
				//triangleCount[std::pair<node,node>(w,u)]++;
				triangleCount[std::pair<node,node>(v,w)]++;
				//triangleCount[std::pair<node,node>(w,v)]++;
			});

			nodeMarker[v] = false;
		});

		removeNode(g, u);
	});

	//The return type of a set is for development purposes only.
	edgeCountSet triangleCountResult;
	triangleCountResult.resize(graph.numberOfEdges());

	for ( auto& kv : triangleCount) {
		triangleCountResult.push_back(std::pair<std::pair<node,node>, count>(kv.first, kv.second));
	}

	return triangleCountResult;
}

void ChibaNishizekiTriangleCounter::removeNode(Graph& graph, const node& u) {
	//isolate the node before removing it.
	graph.forNeighborsOf(u, [&](node v)  {
		graph.removeEdge(u,v);
	});

	graph.removeNode(u);
}

void ChibaNishizekiTriangleCounter::triangleFound(const edgeCountMap& triangleCount, const node& u, const node& v, const node& w) {
	/*triangleCount[std::pair<node,node>(u,v)]++;
	triangleCount[std::pair<node,node>(v,u)]++;
	triangleCount[std::pair<node,node>(u,w)]++;
	triangleCount[std::pair<node,node>(w,u)]++;
	triangleCount[std::pair<node,node>(v,w)]++;
	triangleCount[std::pair<node,node>(w,v)]++;*/
}

} /* namespace NetworKit */
