/*
 * Sparsifiers.cpp
 *
 *  Created on: 24.07.2014
 *      Author: Gerd Lindner
 */

#include "Sparsifiers.h"
#include "../edgeattributes/ChibaNishizekiTriangleCounter.h"
#include "SimmelianJaccardAttributizer.h"
#include "SimmelianOverlapAttributizer.h"
#include "MultiscaleAttributizer.h"
#include "LocalSimilarityAttributizer.h"
#include "RandomEdgeAttributizer.h"
#include "GlobalThresholdFilter.h"

#include "../auxiliary/Random.h"

namespace NetworKit {

/**
 * ---------------------------------------------------------------------------
 */

Sparsifier::Sparsifier(const Graph& inputGraph) : inputGraph(inputGraph) {
}

Graph Sparsifier::getGraph() {
	if (!hasOutput) throw std::runtime_error("Error: run must be called first");

	hasOutput = false;
	return std::move(outputGraph);
}

/**
 * ---------------------------------------------------------------------------
 */

SimmelianBackboneNonParametric::SimmelianBackboneNonParametric(const Graph& graph, double threshold) :
		Sparsifier(graph), threshold(threshold) {}

void SimmelianBackboneNonParametric::run() {
	ChibaNishizekiTriangleCounter triangleAttributizer(inputGraph);
	std::vector<count> triangles = triangleAttributizer.getAttribute();

	SimmelianJaccardAttributizer jaccardAttributizer(inputGraph, triangles);
	std::vector<double> jaccard = jaccardAttributizer.getAttribute();

	GlobalThresholdFilter filter(inputGraph, jaccard, threshold, true);
	outputGraph = filter.calculate();
	hasOutput = true;
}

/**
 * ---------------------------------------------------------------------------
 */

SimmelianBackboneParametric::SimmelianBackboneParametric(const Graph& graph, int maxRank, int minOverlap) :
		Sparsifier(graph), maxRank(maxRank), minOverlap(minOverlap) {}

void SimmelianBackboneParametric::run() {
	ChibaNishizekiTriangleCounter triangleAttributizer(inputGraph);
	std::vector<count> triangles = triangleAttributizer.getAttribute();

	SimmelianOverlapAttributizer overlapAttributizer(inputGraph, triangles, maxRank);
	std::vector<double> overlap = overlapAttributizer.getAttribute();

	GlobalThresholdFilter filter(inputGraph, overlap, minOverlap, true);
	outputGraph = filter.calculate();
	hasOutput = true;
}

/**
 * ---------------------------------------------------------------------------
 */

MultiscaleBackbone::MultiscaleBackbone(const Graph& graph, double alpha) :
		Sparsifier(graph), alpha(alpha) {}

void MultiscaleBackbone::run() {
	std::vector<double> weight(inputGraph.upperEdgeIdBound());
	inputGraph.forEdges([&](node u, node v, edgeid eid) {
		weight[eid] = inputGraph.weight(u, v);
	});

	MultiscaleAttributizer multiscaleAttributizer(inputGraph, weight);
	std::vector<double> multiscale = multiscaleAttributizer.getAttribute();

	GlobalThresholdFilter filter(inputGraph, multiscale, alpha, true);
	outputGraph = filter.calculate();
	hasOutput = true;
}

/**
 * ---------------------------------------------------------------------------
 */

LocalSimilarityBackbone::LocalSimilarityBackbone(const Graph& graph, double e) :
		Sparsifier(graph), e(e) {}

void LocalSimilarityBackbone::run() {
	ChibaNishizekiTriangleCounter triangleCounter(inputGraph);
	std::vector<count> triangles = triangleCounter.getAttribute();

	LocalSimilarityAttributizer localSimAttributizer(inputGraph, triangles);
	std::vector<double> minExponent = localSimAttributizer.getAttribute();

	GlobalThresholdFilter filter(inputGraph, minExponent, e, true);
	outputGraph = filter.calculate();
	hasOutput = true;
}

/**
 * ---------------------------------------------------------------------------
 */

SimmelianMultiscaleBackbone::SimmelianMultiscaleBackbone(const Graph& graph, double alpha) :
		Sparsifier(graph), alpha(alpha) {}

void SimmelianMultiscaleBackbone::run() {
	ChibaNishizekiTriangleCounter triangleAttributizer(inputGraph);
	std::vector<count> triangles = triangleAttributizer.getAttribute();
	std::vector<double> triangles_d = std::vector<double>(triangles.begin(), triangles.end());

	MultiscaleAttributizer multiscaleAttributizer (inputGraph, triangles_d);
	std::vector<double> multiscale = multiscaleAttributizer.getAttribute();

	GlobalThresholdFilter filter(inputGraph, multiscale, alpha, true);
	outputGraph = filter.calculate();
	hasOutput = true;
}

/**
* ---------------------------------------------------------------------------
*/

RandomBackbone::RandomBackbone(const Graph& graph, double ratio) :
		Sparsifier(graph), ratio(ratio) {}

void RandomBackbone::run() {
	RandomEdgeAttributizer randomAttributizer (inputGraph);
	std::vector<double> random = randomAttributizer.getAttribute();

	GlobalThresholdFilter filter(inputGraph, random, ratio, true);
	outputGraph = filter.calculate();
	hasOutput = true;
}

} /* namespace NetworKit */
