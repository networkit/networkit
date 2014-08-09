/*
 * Backbones.cpp
 *
 *  Created on: 24.07.2014
 *      Author: Gerd Lindner
 */

#include "Backbones.h"

namespace NetworKit {

SimmelianBackboneNonParametric::SimmelianBackboneNonParametric(double threshold) : threshold(threshold) {}

Graph SimmelianBackboneNonParametric::calculate(Graph& g, const EdgeAttribute& attribute) {
	g.indexEdges();

	ChibaNishizekiTriangleCounter triangleAttributizer;
	EdgeAttribute triangles = triangleAttributizer.getAttribute(g, EdgeAttribute(g.upperEdgeIdBound()));

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

Graph SimmelianBackboneParametric::calculate(Graph& g, const EdgeAttribute& attribute) {
	g.indexEdges();

	ChibaNishizekiTriangleCounter triangleAttributizer;
	EdgeAttribute triangles = triangleAttributizer.getAttribute(g, EdgeAttribute(g.upperEdgeIdBound()));

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

Graph MultiscaleBackbone::calculate(Graph& g, const EdgeAttribute& attribute) {
	g.indexEdges();

	//TODO: the following will be obsolete once graph edge attributes are used.
	EdgeAttribute weight(g.upperEdgeIdBound());
	g.forEdges([&](node u, node v, edgeid eid) {
		weight[eid] = g.weight(u, v);
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

Graph LocalSimilarityBackbone::calculate(Graph& g, const EdgeAttribute& attribute) {
	g.indexEdges();

	LocalSimilarityAttributizer localSimAttributizer;
	EdgeAttribute minExponent = localSimAttributizer.getAttribute(g, EdgeAttribute(g.upperEdgeIdBound()));

	GlobalThresholdFilter filter(e, true);
	return filter.calculate(g, minExponent);
}

/**
 * ---------------------------------------------------------------------------
 */

SimmelianMultiscaleBackbone::SimmelianMultiscaleBackbone(double alpha) :
		alpha(alpha) {}

Graph SimmelianMultiscaleBackbone::calculate(Graph& g, const EdgeAttribute& attribute) {
	g.indexEdges();

	ChibaNishizekiTriangleCounter triangleAttributizer;
	EdgeAttribute triangles = triangleAttributizer.getAttribute(g, EdgeAttribute(g.upperEdgeIdBound()));

	MultiscaleAttributizer multiscaleAttributizer;
	EdgeAttribute multiscale = multiscaleAttributizer.getAttribute(g, triangles);

	GlobalThresholdFilter filter(alpha, false);
	return filter.calculate(g, multiscale);
}

} /* namespace NetworKit */
