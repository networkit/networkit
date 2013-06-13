/*
 * AlgebraicDistances.h
 *
 *  Created on: Jun 11, 2013
 *      Author: Henning
 */

#ifndef ALGEBRAICDISTANCES_H_
#define ALGEBRAICDISTANCES_H_

#include "../graph/Graph.h"
#include "../auxiliary/RandomProbability.h"
#include <vector>

namespace NetworKit {

const index MAX_NORM = 0;

/**
 * Algebraic distances according to SISC paper by Jie Chen and Ilya Safro
 */
class AlgebraicDistances {
protected:
	std::vector<std::vector<double> > loads; //!< loads[i]: vector of loads of length n for one system
	count numSystems; //!< number of vectors/systems used for algebraic iteration
	count numIters; //!< number of iterations in each system
	const Graph& g; //!< graph on which ADs are computed

	void randomInit();

public:
	AlgebraicDistances(const Graph& graph);
	virtual ~AlgebraicDistances();

	/**
	 * @param numberSystems Number of vectors/systems used for algebraic iteration
	 * @param numberIterations Number of iterations in each system
	 * @param omega Overrelaxation parameter
	 *
	 * Starting with random initialization, compute for all @a numberSystems
	 * "diffusion" systems the situation after @a numberIterations iterations
	 * of overrelaxation with overrelaxation parameter @a omega.
	 *
	 * REQ: Needs to be called before algdist delivers meaningful results!
	 */
	void preprocess(count numberSystems, count numberIterations, double omega);

	/**
	 * @return Extended algebraic distance between node @a u and node @a v in norm @a norm with
	 * default norm 2.
	 *
	 * Maximum norm is realized by setting @a norm to 0.
	 */
	double algdist(node u, node v, index norm = 2) const;

	/**
	 * @return Geometric mean of "load" values of node @a u
	 */
	double geometricMeanLoad(node u) const;
};

} /* namespace NetworKit */

#endif /* ALGEBRAICDISTANCES_H_ */
