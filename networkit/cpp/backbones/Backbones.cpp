/*
 * Backbones.cpp
 *
 *  Created on: 24.07.2014
 *      Author: Gerd Lindner
 */

#include "Backbones.h"
#include "../edgeproperties/ChibaNishizekiTriangleCounter.h"
#include "SimmelianJaccardAttributizer.h"
#include "SimmelianOverlapAttributizer.h"
#include "MultiscaleAttributizer.h"
#include "LocalSimilarityAttributizer.h"
#include "RandomAttributizer.h"
#include "GlobalThresholdFilter.h"

#include "../auxiliary/Random.h"

namespace NetworKit {

SimmelianBackboneNonParametric::SimmelianBackboneNonParametric(double threshold) : threshold(threshold) {}

Graph SimmelianBackboneNonParametric::calculate(const Graph& g) {
	ChibaNishizekiTriangleCounter triangleAttributizer(g);
	std::vector<count> triangles = triangleAttributizer.getAttribute();

	SimmelianJaccardAttributizer jaccardAttributizer(g, triangles);
	std::vector<double> jaccard = jaccardAttributizer.getAttribute();

	GlobalThresholdFilter filter(g, jaccard, threshold, true);
	return filter.calculate();
}

/**
 * ---------------------------------------------------------------------------
 */

SimmelianBackboneParametric::SimmelianBackboneParametric(int maxRank, int minOverlap) :
		maxRank(maxRank), minOverlap(minOverlap) {}

Graph SimmelianBackboneParametric::calculate(const Graph& g) {
	ChibaNishizekiTriangleCounter triangleAttributizer(g);
	std::vector<count> triangles = triangleAttributizer.getAttribute();

	SimmelianOverlapAttributizer overlapAttributizer(g, triangles, maxRank);
	std::vector<double> overlap = overlapAttributizer.getAttribute();

	GlobalThresholdFilter filter(g, overlap, minOverlap, true);
	return filter.calculate();
}

/**
 * ---------------------------------------------------------------------------
 */

MultiscaleBackbone::MultiscaleBackbone(double alpha) :
		alpha(alpha) {}

Graph MultiscaleBackbone::calculate(const Graph& g) {
	std::vector<double> weight(g.upperEdgeIdBound());
	g.forEdges([&](node u, node v, edgeid eid) {
		weight[eid] = g.weight(u, v);
	});

	MultiscaleAttributizer multiscaleAttributizer(g, weight);
	std::vector<double> multiscale = multiscaleAttributizer.getAttribute();

	GlobalThresholdFilter filter(g, multiscale, alpha, false);
	return filter.calculate();
}

/**
 * ---------------------------------------------------------------------------
 */

LocalSimilarityBackbone::LocalSimilarityBackbone(double e) :
		e(e) {}

Graph LocalSimilarityBackbone::calculate(const Graph& g) {
	ChibaNishizekiTriangleCounter triangleCounter(g);
	std::vector<count> triangles = triangleCounter.getAttribute();

	LocalSimilarityAttributizer localSimAttributizer(g, triangles);
	std::vector<double> minExponent = localSimAttributizer.getAttribute();

	GlobalThresholdFilter filter(g, minExponent, e, false);
	return filter.calculate();
}

/**
 * ---------------------------------------------------------------------------
 */

SimmelianMultiscaleBackbone::SimmelianMultiscaleBackbone(double alpha) :
		alpha(alpha) {}

Graph SimmelianMultiscaleBackbone::calculate(const Graph& g) {
	ChibaNishizekiTriangleCounter triangleAttributizer(g);
	std::vector<count> triangles = triangleAttributizer.getAttribute();
	std::vector<double> triangles_d = std::vector<double>(triangles.begin(), triangles.end());

	MultiscaleAttributizer multiscaleAttributizer (g, triangles_d);
	std::vector<double> multiscale = multiscaleAttributizer.getAttribute();

	GlobalThresholdFilter filter(g, multiscale, alpha, false);
	return filter.calculate();
}

/**
* ---------------------------------------------------------------------------
*/

RandomBackbone::RandomBackbone(double ratio) :
		ratio(ratio) {}

Graph RandomBackbone::calculate(const Graph& g) {
	RandomAttributizer randomAttributizer (g);
	std::vector<double> random = randomAttributizer.getAttribute();

	GlobalThresholdFilter filter(g, random, ratio, false);
	return filter.calculate();
}

} /* namespace NetworKit */
