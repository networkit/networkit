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
#include "ChibaNishizekiTriangleCounter.h"
#include "SimmelianJaccardAttributizer.h"
#include "SimmelianOverlapAttributizer.h"
#include "BackboneCalculator.h"
#include "GlobalThresholdFilter.h"

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

class SimmelianBackboneNonparametric : public BackboneCalculator {

public:
	/**
	 * Creates a new instance of the non-parametric (jaccard) variant of the Simmelian Backbone calculator.
	 * @param threshold		the jaccard index threshold.
	 */
	SimmelianBackboneNonparametric(double threshold);

	Graph calculate(const Graph& graph, const edgeAttribute& attribute);

private:
	double threshold;

};

SimmelianBackboneNonparametric::SimmelianBackboneNonparametric(double threshold) : threshold(threshold) {}

Graph SimmelianBackboneNonparametric::calculate(const Graph& g, const edgeAttribute& attribute) {
	ChibaNishizekiTriangleCounter triangleAttributizer;
	edgeAttribute triangles = triangleAttributizer.getAttribute(g, edgeAttribute(0));

	SimmelianJaccardAttributizer jaccardAttributizer;
	edgeAttribute jaccard = jaccardAttributizer.getAttribute(g, triangles);

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

	Graph calculate(const Graph& graph, const edgeAttribute& attribute);

private:
	int maxRank;
	int minOverlap;

};

SimmelianBackboneParametric::SimmelianBackboneParametric(int maxRank, int minOverlap) :
		maxRank(maxRank), minOverlap(minOverlap) {}

Graph SimmelianBackboneParametric::calculate(const Graph& g, const edgeAttribute& attribute) {
	ChibaNishizekiTriangleCounter triangleAttributizer;
	edgeAttribute triangles = triangleAttributizer.getAttribute(g, edgeAttribute(0));

	SimmelianOverlapAttributizer overlapAttributizer(maxRank);
	edgeAttribute overlap = overlapAttributizer.getAttribute(g, triangles);

	GlobalThresholdFilter filter(minOverlap, true);
	return filter.calculate(g, overlap);
}


} /* namespace NetworKit */
#endif /* BACKBONECALCULATOR_H_ */
