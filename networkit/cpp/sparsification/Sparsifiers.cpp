/*
 * Sparsifiers.cpp
 *
 *  Created on: 24.07.2014
 *      Author: Gerd Lindner
 */

#include "Sparsifiers.h"
#include "../edgescores/ChibaNishizekiTriangleEdgeScore.h"
#include "../edgescores/PrefixJaccardScore.h"
#include "SimmelianOverlapScore.h"
#include "MultiscaleScore.h"
#include "LocalSimilarityScore.h"
#include "RandomEdgeScore.h"
#include "GlobalThresholdFilter.h"

#include "../auxiliary/Random.h"

namespace NetworKit {

Sparsifier::Sparsifier(const Graph& inputGraph) : inputGraph(inputGraph) {
}

Graph Sparsifier::getGraph() {
	if (!hasOutput) throw std::runtime_error("Error: run must be called first");

	hasOutput = false;
	return std::move(outputGraph);
}


SimmelianBackboneNonParametric::SimmelianBackboneNonParametric(const Graph& graph, double threshold) :
		Sparsifier(graph), threshold(threshold) {}

void SimmelianBackboneNonParametric::run() {
	ChibaNishizekiTriangleEdgeScore triangleEdgeScore(inputGraph);
	triangleEdgeScore.run();
	std::vector<count> triangles = triangleEdgeScore.scores();

	PrefixJaccardScore<count> jaccardScore(inputGraph, triangles);
	jaccardScore.run();
	std::vector<double> jaccard = jaccardScore.scores();

	GlobalThresholdFilter filter(inputGraph, jaccard, threshold, true);
	outputGraph = filter.calculate();
	hasOutput = true;
}


SimmelianBackboneParametric::SimmelianBackboneParametric(const Graph& graph, int maxRank, int minOverlap) :
		Sparsifier(graph), maxRank(maxRank), minOverlap(minOverlap) {}

void SimmelianBackboneParametric::run() {
	ChibaNishizekiTriangleEdgeScore triangleEdgeScore(inputGraph);
	triangleEdgeScore.run();
	std::vector<count> triangles = triangleEdgeScore.scores();

	SimmelianOverlapScore overlapScore(inputGraph, triangles, maxRank);
	overlapScore.run();
	std::vector<double> overlap = overlapScore.scores();

	GlobalThresholdFilter filter(inputGraph, overlap, minOverlap, true);
	outputGraph = filter.calculate();
	hasOutput = true;
}


MultiscaleBackbone::MultiscaleBackbone(const Graph& graph, double alpha) :
		Sparsifier(graph), alpha(alpha) {}

void MultiscaleBackbone::run() {
	std::vector<double> weight(inputGraph.upperEdgeIdBound());
	inputGraph.forEdges([&](node u, node v, edgeid eid) {
		weight[eid] = inputGraph.weight(u, v);
	});

	MultiscaleScore multiscaleScorer(inputGraph, weight);
	multiscaleScorer.run();
	std::vector<double> multiscale = multiscaleScorer.scores();

	GlobalThresholdFilter filter(inputGraph, multiscale, alpha, true);
	outputGraph = filter.calculate();
	hasOutput = true;
}


LocalSimilarityBackbone::LocalSimilarityBackbone(const Graph& graph, double e) :
		Sparsifier(graph), e(e) {}

void LocalSimilarityBackbone::run() {
	ChibaNishizekiTriangleEdgeScore triangleEdgeScore(inputGraph);
	triangleEdgeScore.run();
	std::vector<count> triangles = triangleEdgeScore.scores();

	LocalSimilarityScore localSimScore(inputGraph, triangles);
	localSimScore.run();
	std::vector<double> minExponent = localSimScore.scores();

	GlobalThresholdFilter filter(inputGraph, minExponent, e, true);
	outputGraph = filter.calculate();
	hasOutput = true;
}

SimmelianMultiscaleBackbone::SimmelianMultiscaleBackbone(const Graph& graph, double alpha) :
		Sparsifier(graph), alpha(alpha) {}

void SimmelianMultiscaleBackbone::run() {
	ChibaNishizekiTriangleEdgeScore triangleEdgeScore(inputGraph);
	triangleEdgeScore.run();
	std::vector<count> triangles = triangleEdgeScore.scores();
	std::vector<double> triangles_d = std::vector<double>(triangles.begin(), triangles.end());

	MultiscaleScore multiscaleScorer (inputGraph, triangles_d);
	multiscaleScorer.run();
	std::vector<double> multiscale = multiscaleScorer.scores();

	GlobalThresholdFilter filter(inputGraph, multiscale, alpha, true);
	outputGraph = filter.calculate();
	hasOutput = true;
}

RandomBackbone::RandomBackbone(const Graph& graph, double ratio) :
		Sparsifier(graph), ratio(ratio) {}

void RandomBackbone::run() {
	RandomEdgeScore randomScorer (inputGraph);
	randomScorer.run();
	std::vector<double> random = randomScorer.scores();

	GlobalThresholdFilter filter(inputGraph, random, ratio, true);
	outputGraph = filter.calculate();
	hasOutput = true;
}

} /* namespace NetworKit */
