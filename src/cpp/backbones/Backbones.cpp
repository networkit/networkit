/*
 * Backbones.cpp
 *
 *  Created on: 24.07.2014
 *      Author: Gerd Lindner
 */

#include "Backbones.h"

namespace NetworKit {

SimmelianBackboneNonParametric::SimmelianBackboneNonParametric(double threshold) : threshold(threshold) {}

Graph SimmelianBackboneNonParametric::calculate(const Graph& g, const EdgeAttribute& attribute) {
	ChibaNishizekiTriangleCounter triangleAttributizer;
	EdgeAttribute triangles = triangleAttributizer.getAttribute(g, EdgeAttribute());

	SimmelianJaccardAttributizer jaccardAttributizer;
	EdgeAttribute jaccard = jaccardAttributizer.getAttribute(g, triangles);

	GlobalThresholdFilter filter(threshold, true);
	return filter.calculate(g, jaccard);
}

/**
 * ---------------------------------------------------------------------------
 */

SimmelianBackboneParametric::SimmelianBackboneParametric(int maxRank, int minOverlap) :
		maxRank(maxRank), minOverlap(minOverlap) {}

Graph SimmelianBackboneParametric::calculate(const Graph& g, const EdgeAttribute& attribute) {
	ChibaNishizekiTriangleCounter triangleAttributizer;
	EdgeAttribute triangles = triangleAttributizer.getAttribute(g, EdgeAttribute());

	SimmelianOverlapAttributizer overlapAttributizer(maxRank);
	EdgeAttribute overlap = overlapAttributizer.getAttribute(g, triangles);

	GlobalThresholdFilter filter(minOverlap, true);
	return filter.calculate(g, overlap);
}

/**
 * ---------------------------------------------------------------------------
 */

MultiscaleBackbone::MultiscaleBackbone(double alpha) :
		alpha(alpha) {}

Graph MultiscaleBackbone::calculate(const Graph& g, const EdgeAttribute& attribute) {
	//TODO: the following will be obsolete once graph edge attributes are used.
	EdgeAttribute weight;
	g.forEdges([&](node u, node v) {
		weight.set(uEdge(u, v), g.weight(u, v));
	});

	MultiscaleAttributizer multiscaleAttributizer;
	EdgeAttribute multiscale = multiscaleAttributizer.getAttribute(g, weight);

	GlobalThresholdFilter filter(alpha, false);
	return filter.calculate(g, multiscale);
}

/**
 * ---------------------------------------------------------------------------
 */

LocalSimilarityBackbone::LocalSimilarityBackbone(double e) :
		e(e) {}

Graph LocalSimilarityBackbone::calculate(const Graph& g, const EdgeAttribute& attribute) {
	LocalSimilarityAttributizer localSimAttributizer;
	EdgeAttribute minExponent = localSimAttributizer.getAttribute(g, EdgeAttribute());

	GlobalThresholdFilter filter(e, true);
	return filter.calculate(g, minExponent);
}

/**
 * ---------------------------------------------------------------------------
 */

SimmelianMultiscaleBackbone::SimmelianMultiscaleBackbone(double alpha) :
		alpha(alpha) {}

Graph SimmelianMultiscaleBackbone::calculate(const Graph& g, const EdgeAttribute& attribute) {
	ChibaNishizekiTriangleCounter triangleAttributizer;
	EdgeAttribute triangles = triangleAttributizer.getAttribute(g, EdgeAttribute());

	MultiscaleAttributizer multiscaleAttributizer;
	EdgeAttribute multiscale = multiscaleAttributizer.getAttribute(g, triangles);

	GlobalThresholdFilter filter(alpha, false);
	return filter.calculate(g, multiscale);
}

} /* namespace NetworKit */
