/*
 * Sparsifiers.h
 *
 *  Created on: 23.07.2014
 *      Author: Gerd Lindner
 */

#ifndef SPARSIFIERS_H_
#define SPARSIFIERS_H_

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Here, we combine attribute generators and edge attribute filters into different
 * backbone algorithms.
 */

/**
 * Abstract base class for Sparsifiers.
 */
class Sparsifier {

public:

	Sparsifier(const Graph& inputGraph);

	virtual ~Sparsifier() = default;

	/**
	 * REQ: Needs to fill outputGraph and set hasOutput to true.
	 */
	virtual void run() = 0;

	Graph getGraph();

protected:
	const Graph& inputGraph;
	Graph outputGraph;
	bool hasOutput;

};


/**
 * --------------------------------------------------------------------------------------
 * Simmelian Backbone: Non-parametric variant (Jaccard)
 * --------------------------------------------------------------------------------------
 */

class SimmelianBackboneNonParametric : public Sparsifier {

public:
	/**
	 * Creates a new instance of the non-parametric (jaccard) variant of the Simmelian Backbone calculator.
	 * @param graph			the input graph
	 * @param threshold		the jaccard index threshold.
	 */
	SimmelianBackboneNonParametric(const Graph& graph, double threshold);

	virtual void run() override;

private:
	double threshold;

};

/**
 * --------------------------------------------------------------------------------------
 * Simmelian Backbone: Parametric variant (Top-k neighborhood overlap)
 * --------------------------------------------------------------------------------------
 */

class SimmelianBackboneParametric : public Sparsifier {

public:
	/**
	 * Creates a new instance of the parametric variant of the Simmelian Backbone calculator.
	 * @param graph			the input graph
	 * @param maxRank 		the maximum rank that is considered for overlap calculation
	 * @param minOverlap	the minimum overlap of the top-k neighbors for an edge to be
		 	 	 	 	 	 	contained in the backbone.
	 */
	SimmelianBackboneParametric(const Graph& graph, int maxRank, int minOverlap);

	virtual void run() override;

private:
	int maxRank;
	int minOverlap;

};

/**
 * --------------------------------------------------------------------------------------
 * Multiscale Backbone
 * --------------------------------------------------------------------------------------
 */

class MultiscaleBackbone : public Sparsifier {

public:
	/**
	 * Creates a new instance of the parametric variant of the Simmelian Backbone calculator.
	 * @param graph			the input graph
	 * @param alpha 		the probability threshold
	 */
	MultiscaleBackbone(const Graph& graph, double alpha);

	virtual void run() override;

private:
	double alpha;

};

/**
 * --------------------------------------------------------------------------------------
 * Local Similarity Backbone
 * --------------------------------------------------------------------------------------
 */

class LocalSimilarityBackbone : public Sparsifier {

public:
	/**
	 * Creates a new instance of the Local Similarity Backbone calculator
	 * @param graph			the input graph
	 * @param e				the threshold value
	 */
	LocalSimilarityBackbone(const Graph& graph, double e);

	virtual void run() override;

private:
	double e;

};

/**
 * --------------------------------------------------------------------------------------
 * Multiscale backbone using simmelianness as weight
 * --------------------------------------------------------------------------------------
 */

class SimmelianMultiscaleBackbone : public Sparsifier {

public:
	/**
	 * Creates a new instance of the Simmelian Multiscale Backbone calculator
	 * @param graph			the input graph
	 * @param alpha			the threshold value for multiscale filtering
	 */
	SimmelianMultiscaleBackbone(const Graph& graph, double alpha);

	virtual void run() override;

private:
	double alpha;

};

/**
* --------------------------------------------------------------------------------------
* Backbone that contains approximately a given percentage
* of edges of the original graph (mainly for comparison purposes).
* The edges are selected unformly at random.
* --------------------------------------------------------------------------------------
*/

class RandomBackbone : public Sparsifier {

public:
	/**
	* Creates a new instance of the Random Backbone calculator
	* @param graph			the input graph
	* @param ratio			edge ratio in [0,1] to be kept in backbone.
	*/
	RandomBackbone(const Graph& graph, double ratio);

	virtual void run() override;

private:
	double ratio;

};



} /* namespace NetworKit */
#endif /* SPARSIFIERS_H_ */
