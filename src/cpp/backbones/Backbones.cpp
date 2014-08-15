/*
 * Backbones.cpp
 *
 *  Created on: 24.07.2014
 *      Author: Gerd Lindner
 */

#include "Backbones.h"
#include "ChibaNishizekiTriangleCounter.h"
#include "SimmelianJaccardAttributizer.h"
#include "SimmelianOverlapAttributizer.h"
#include "MultiscaleAttributizer.h"
#include "LocalSimilarityAttributizer.h"
#include "RandomAttributizer.h"
#include "GlobalThresholdFilter.h"

#include "../auxiliary/Random.h"

namespace NetworKit {

SimmelianBackboneNonParametric::SimmelianBackboneNonParametric(double threshold) : threshold(threshold) {}

Graph SimmelianBackboneNonParametric::calculate(Graph& g) {
	g.indexEdges();

	ChibaNishizekiTriangleCounter triangleAttributizer;
	std::vector<int> triangles = triangleAttributizer.getAttribute(g, std::vector<int>(g.upperEdgeIdBound()));

	SimmelianJaccardAttributizer jaccardAttributizer;
	std::vector<double> jaccard = jaccardAttributizer.getAttribute(g, triangles);

	GlobalThresholdFilter filter(threshold, true);
	return filter.calculate(g, jaccard);
}

/**
 * ---------------------------------------------------------------------------
 */

SimmelianBackboneParametric::SimmelianBackboneParametric(int maxRank, int minOverlap) :
		maxRank(maxRank), minOverlap(minOverlap) {}

Graph SimmelianBackboneParametric::calculate(Graph& g) {
	g.indexEdges();

	ChibaNishizekiTriangleCounter triangleAttributizer;
	std::vector<int> triangles = triangleAttributizer.getAttribute(g, std::vector<int>(g.upperEdgeIdBound()));

	SimmelianOverlapAttributizer overlapAttributizer(maxRank);
	std::vector<double> overlap = overlapAttributizer.getAttribute(g, triangles);

	GlobalThresholdFilter filter(minOverlap, true);
	return filter.calculate(g, overlap);
}

/**
 * ---------------------------------------------------------------------------
 */

MultiscaleBackbone::MultiscaleBackbone(double alpha) :
		alpha(alpha) {}

Graph MultiscaleBackbone::calculate(Graph& g) {
	g.indexEdges();

	//TODO: the following will be obsolete once graph edge attributes are used. ?
	std::vector<double> weight(g.upperEdgeIdBound());
	g.forEdges([&](node u, node v, edgeid eid) {
		weight[eid] = g.weight(u, v);
	});

	MultiscaleAttributizer multiscaleAttributizer;
	std::vector<double> multiscale = multiscaleAttributizer.getAttribute(g, weight);

	GlobalThresholdFilter filter(alpha, false);
	return filter.calculate(g, multiscale);
}

/**
 * ---------------------------------------------------------------------------
 */

LocalSimilarityBackbone::LocalSimilarityBackbone(double e) :
		e(e) {}

Graph LocalSimilarityBackbone::calculate(Graph& g) {
	g.indexEdges();

	LocalSimilarityAttributizer localSimAttributizer;
	std::vector<double> minExponent = localSimAttributizer.getAttribute(g, std::vector<int>(g.upperEdgeIdBound()));

	GlobalThresholdFilter filter(e, true);
	return filter.calculate(g, minExponent);
}

/**
 * ---------------------------------------------------------------------------
 */

SimmelianMultiscaleBackbone::SimmelianMultiscaleBackbone(double alpha) :
		alpha(alpha) {}

Graph SimmelianMultiscaleBackbone::calculate(Graph& g) {
	g.indexEdges();

	ChibaNishizekiTriangleCounter triangleAttributizer;
	std::vector<int> triangles = triangleAttributizer.getAttribute(g, std::vector<int>(g.upperEdgeIdBound()));

	MultiscaleAttributizer multiscaleAttributizer;
	std::vector<double> multiscale = multiscaleAttributizer.getAttribute(g,
			std::vector<double>(triangles.begin(), triangles.end()));

	GlobalThresholdFilter filter(alpha, false);
	return filter.calculate(g, multiscale);
}

/**
* ---------------------------------------------------------------------------
*/

RandomBackbone::RandomBackbone(double ratio) :
		ratio(ratio) {}

Graph RandomBackbone::calculate(Graph& g) {
	g.indexEdges();

	RandomAttributizer randomAttributizer;
	std::vector<double> random = randomAttributizer.getAttribute(g, std::vector<int>(g.upperEdgeIdBound()));

	GlobalThresholdFilter filter(ratio, false);
	return filter.calculate(g, random);
} /* namespace NetworKit */

}
