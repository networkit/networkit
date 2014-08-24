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

std::vector<int> ChibaNishizekiTriangleCounter::getAttribute(const Graph& graph, const std::vector<int>& attribute) {
	Graph g = graph;

	//Node attribute: marker
	std::vector<bool> nodeMarker(graph.upperNodeIdBound(), false);

	//Edge attribute: triangle count
	std::vector<int> triangleCount(graph.upperEdgeIdBound(), 0);

	g.forNodes([&](node u) {
		//Mark all neighbors
		g.forNeighborsOf(u, [&](node v) {
			nodeMarker[v] = true;
		});

		//TODO: Remove _u currently breaks everything due to missing lambda case in Graph.h
		//For all neighbors: check for already marked neighbors.
		g.forNeighborsOf(u, [&](node _u, node v, edgeid eid_uv) {
			g.forNeighborsOf(v, [&](node _v, node w, edgeid eid_vw) {
				if (nodeMarker[w]) {

					edgeid eid_uw = graph.edgeId(u, w);

					triangleCount[eid_uv] = triangleCount[eid_uv] + 1;
					triangleCount[eid_uw] = triangleCount[eid_uw] + 1;
					triangleCount[eid_vw] = triangleCount[eid_vw] + 1;
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

} /* namespace NetworKit */
