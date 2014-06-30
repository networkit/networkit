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

edgeCountMap ChibaNishizekiTriangleCounter::triangleCounts(const Graph& graph) {
	INFO("Started triangle counting...", "\n");

	Aux::Timer timer;
	timer.start();

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

	timer.stop();
	INFO("elapsed millisecs for CN triangle counting: ", timer.elapsedMilliseconds(), "\n");

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

} /* namespace NetworKit */
