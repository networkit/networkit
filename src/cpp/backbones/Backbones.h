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

	Graph calculate(const Graph& graph);

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

	Graph calculate(const Graph& graph);

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

	Graph calculate(const Graph& graph);

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

	Graph calculate(const Graph& graph);

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

	Graph calculate(const Graph& graph);

private:
	double alpha;

};

/**
* --------------------------------------------------------------------------------------
* Backbone that contains approximately a given percentage
* of edges of the original graph (mainly for comparison purposes).
* The edges are selected randomly.
* --------------------------------------------------------------------------------------
*/

class RandomBackbone : public BackboneCalculator {

public:
	/**
		* Creates a new instance of the Random Backbone calculator
		* @param ratio			edge ratio in [0,1] to be kept in backbone.
		*/
	RandomBackbone(double ratio);

	Graph calculate(const Graph& graph);

private:
	double ratio;

};



} /* namespace NetworKit */
#endif /* BACKBONECALCULATOR_H_ */
