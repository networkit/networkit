/*
 * AlgebraicDistance.h
 *
 *  Created on: 19.06.2013
 *      Author: cls
 */

#ifndef ALGEBRAICDISTANCE_H_
#define ALGEBRAICDISTANCE_H_

#include "NodeDistance.h"
#include "../graph/Graph.h"
#include "../auxiliary/RandomProbability.h"

namespace NetworKit {


class AlgebraicDistance: public NetworKit::NodeDistance {
public:

	/**
	 * @param numberSystems Number of vectors/systems used for algebraic iteration
	 * @param numberIterations Number of iterations in each system
	 * @param omega Overrelaxation parameter
	 */
	AlgebraicDistance(const Graph& G, count numberSystems, count numberIterations, double omega = 0.5, index norm = 2);

	 ~AlgebraicDistance();

	/**
	 * Starting with random initialization, compute for all @a numberSystems
	 * "diffusion" systems the situation after @a numberIterations iterations
	 * of overrelaxation with overrelaxation parameter @a omega.
	 *
	 * REQ: Needs to be called before algdist delivers meaningful results!
	 */
	 virtual void preprocess();

	/**
	 * @return Extended algebraic distance between node @a u and node @a v in norm @a norm with
	 * default norm 2.
	 *
	 * Maximum norm is realized by setting @a norm to 0.
	 */
	 virtual double distance(node u, node v);

protected:

	count numSystems; //!< number of vectors/systems used for algebraic iteration
	count numIters; //!< number of iterations in each system
	double omega; //!<
	index norm;
	const index MAX_NORM = 0;

	std::vector<std::vector<double> > loads; //!< loads[i]: vector of loads of length n for one system


	void randomInit();
};

} /* namespace NetworKit */
#endif /* ALGEBRAICDISTANCE_H_ */
