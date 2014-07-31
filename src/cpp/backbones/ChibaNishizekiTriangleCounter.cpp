/*
 * ChibaNishizekiTriangleCounter.cpp
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner
 */

#include "ChibaNishizekiTriangleCounter.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Timer.h"

namespace NetworKit {

EdgeAttribute ChibaNishizekiTriangleCounter::getAttribute(const Graph& graph, const EdgeAttribute& attribute) {
	Graph g = graph;

	//Node attribute: marker
	std::vector<bool> nodeMarker(false);
	nodeMarker.resize(graph.numberOfNodes());

	//Edge attribute: triangle count
	EdgeAttribute triangleCount(graph.upperedgeIdBound(), 0.0);

	g.forNodes([&](node u) {
		//Mark all neighbors
		g.forNeighborsOf(u, [&](node v) {
			nodeMarker[v] = true;
		});

		//For all neighbors: check for already marked neighbors.
		g.forNeighborsOf(u, [&](node v) {
			g.forNeighborsOf(v, [&](node w) {
				if (nodeMarker[w]) {
					triangleFound(graph, triangleCount, u, v, w);
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
		const Graph& graph, EdgeAttribute& triangleCount,
		const node& u, const node& v,const node& w ) {
	triangleCount[graph.edgeId(u, v)] = triangleCount[graph.edgeId(u, v)] + 1.0);
	triangleCount[graph.edgeId(u, w)] = triangleCount[graph.edgeId(u, w)] + 1.0);
	triangleCount[graph.edgeId(v, w)] = triangleCount[graph.edgeId(v, w)] + 1.0);
}

} /* namespace NetworKit */
