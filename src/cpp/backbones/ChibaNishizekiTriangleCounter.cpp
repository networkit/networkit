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
	std::vector<bool> nodeMarker(graph.numberOfNodes(), false);

	//Edge attribute: triangle count
	EdgeAttribute triangleCount(graph.upperEdgeIdBound(), 0.0);

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

void ChibaNishizekiTriangleCounter::removeNode(Graph& graph, node u) {
	//isolate the node before removing it.
	graph.forNeighborsOf(u, [&](node v) {
		graph.removeEdge(u,v);
	});

	graph.removeNode(u);
}

void ChibaNishizekiTriangleCounter::triangleFound(
		const Graph& graph, EdgeAttribute& triangleCount,
		node u, node v, node w ) {
	edgeid id_uv = graph.edgeId(u, v);
	edgeid id_uw = graph.edgeId(u, w);
	edgeid id_vw = graph.edgeId(v, w);
	triangleCount[id_uv] = triangleCount[id_uv] + 1.0;
	triangleCount[id_uw] = triangleCount[id_uw] + 1.0;
	triangleCount[id_vw] = triangleCount[id_vw] + 1.0;
}

} /* namespace NetworKit */
