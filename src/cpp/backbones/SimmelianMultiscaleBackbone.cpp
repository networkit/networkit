/*
 * SimmelianMultiscaleBackbone.cpp
 *
 *  Created on: 17.07.2014
 *      Author: Gerd Lindner
 */

#include "SimmelianMultiscaleBackbone.h"
#include "ChibaNishizekiTriangleCounter.h"
#include "MultiscaleBackbone.h"

namespace NetworKit {

SimmelianMultiscaleBackbone::SimmelianMultiscaleBackbone(double threshold) : threshold(threshold) {}

Graph SimmelianMultiscaleBackbone::calculate(const Graph& graph) {
	/**
	 * We create a weighted graph by using triangle counts as edge weights
	 * and then apply the multiscale backbone algorithm.
	 */
	ChibaNishizekiTriangleCounter counter;
	edgeCountMap triangles = counter.triangleCounts(graph);

	//TODO: avoid copying the graph...
	Graph graphCopy = cloneNodes(graph, true);

	graph.forEdges([&](node u, node v) {
		graphCopy.setWeight(u, v, triangles[uEdge(u, v)]);
	});

	MultiscaleBackbone multiscale(threshold);
	return multiscale.calculate(graphCopy);
}

} /* namespace NetworKit */
