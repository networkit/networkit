/*
 * Backbones.h
 *
 *  Created on: 23.07.2014
 *      Author: Gerd Lindner
 */

#ifndef BACKBONES_H_
#define BACKBONES_H_

#include "../graph/Graph.h"
#include "AttributeGenerator.h"
#include "BackboneCalculator.h"
#include "ChibaNishizekiTriangleCounter.h"
#include "SimmelianJaccardAttributizer.h"
#include "SimmelianOverlapAttributizer.h"
#include "GlobalThresholdFilter.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

/** 
 * Here, we combine attribute generators and edge attribute filters into different
 * backbone algorithms.
 */


/**
 * --------------------------------------------------------------------------------------
 * Simmelian Backbone: Non-parametric variant (Jaccard)
 * --------------------------------------------------------------------------------------
 */

class SimmelianBackboneNonParametric : public BackboneCalculator {

public:
	/**
	 * Creates a new instance of the non-parametric (jaccard) variant of the Simmelian Backbone calculator.
	 * @param threshold		the jaccard index threshold.
	 */
	SimmelianBackboneNonParametric(double threshold);

	Graph calculate(const Graph& graph, const EdgeAttribute& attribute);

private:
	double threshold;

};

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
 * --------------------------------------------------------------------------------------
 * Simmelian Backbone: Parametric variant (Top-k neighborhood overlap)
 * --------------------------------------------------------------------------------------
 */

class SimmelianBackboneParametric : public BackboneCalculator {

public:
	/**
		 * Creates a new instance of the parametric variant of the Simmelian Backbone calculator.
		 * @param maxRank 		the maximum rank that is considered for overlap calculation
		 * @param minOverlap	the minimum overlap of the top-k neighbors for an edge to be
		 	 	 	 	 	 	 contained in the backbone.
		 */
	SimmelianBackboneParametric(int maxRank, int minOverlap);

	Graph calculate(const Graph& graph, const EdgeAttribute& attribute);

private:
	int maxRank;
	int minOverlap;

};

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
 * --------------------------------------------------------------------------------------
 * Multiscale Backbone
 * --------------------------------------------------------------------------------------
 */

class MultiscaleBackbone : public BackboneCalculator {

public:
	/**
		 * Creates a new instance of the parametric variant of the Simmelian Backbone calculator.
		 * @param alpha 		the probability threshold
		 */
	MultiscaleBackbone(double alpha);

	Graph calculate(const Graph& graph, const EdgeAttribute& attribute);

private:
	double alpha;

};

MultiscaleBackbone::MultiscaleBackbone(double alpha) :
		alpha(alpha) {}

Graph MultiscaleBackbone::calculate(const Graph& g, const EdgeAttribute& attribute) {
	//TODO: the following will be obsolete once graph edge attributes are used.
	EdgeAttribute weight;
	graph.forEdges([&](node u, node v) {
		weight.set(uEdge(u, v), g.weight(u, v));
	});

	MultiscaleAttributizer multiscaleAttributizer;
	EdgeAttribute multiscale = multiscaleAttributizer.getAttribute(weight);

	GlobalThresholdFilter filter(alpha, false);
	return filter.calculate(g, multiscale);
}



} /* namespace NetworKit */
#endif /* BACKBONECALCULATOR_H_ */
