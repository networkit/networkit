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
#include "MultiscaleAttributizer.h"
#include "LocalSimilarityAttributizer.h"
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

	Graph calculate(Graph& graph, const EdgeAttribute& attribute);

private:
	double threshold;

};

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

	Graph calculate(Graph& graph, const EdgeAttribute& attribute);

private:
	int maxRank;
	int minOverlap;

};

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

	Graph calculate(Graph& graph, const EdgeAttribute& attribute);

private:
	double alpha;

};

/**
 * --------------------------------------------------------------------------------------
 * Local Similarity Backbone
 * --------------------------------------------------------------------------------------
 */

class LocalSimilarityBackbone : public BackboneCalculator {

public:
	/**
		 * Creates a new instance of the Local Similarity Backbone calculator
		 * @param e			the threshold value
		 */
	LocalSimilarityBackbone(double e);

	Graph calculate(Graph& graph, const EdgeAttribute& attribute);

private:
	double e;

};

/**
 * --------------------------------------------------------------------------------------
 * Multiscale backbone using simmelianness as weight
 * --------------------------------------------------------------------------------------
 */

class SimmelianMultiscaleBackbone : public BackboneCalculator {

public:
	/**
		 * Creates a new instance of the Simmelian Multiscale Backbone calculator
		 * @param alpha			the threshold value for multiscale filtering
		 */
	SimmelianMultiscaleBackbone(double alpha);

	Graph calculate(Graph& graph, const EdgeAttribute& attribute);

private:
	double alpha;

};



} /* namespace NetworKit */
#endif /* BACKBONECALCULATOR_H_ */
