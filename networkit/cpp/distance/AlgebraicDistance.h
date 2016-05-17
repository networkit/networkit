/*
 * AlgebraicDistance.h
 *
 *  Created on: 03.11.2015
 *      Author: Henning Meyerhenke, Christian Staudt, Michael Hamann
 */

#ifndef ALGEBRAICDISTANCE_H_
#define ALGEBRAICDISTANCE_H_

#include "NodeDistance.h"


namespace NetworKit {

/**
 * @ingroup distance
 * Algebraic distance assigns a distance value to pairs of nodes
 * according to their structural closeness in the graph.
 * Algebraic distances will become small within dense subgraphs.
 *
 */
class AlgebraicDistance: public NetworKit::NodeDistance {

public:

	/**
	 * @param G The graph.
	 * @param numberSystems Number of vectors/systems used for algebraic iteration.
	 * @param numberIterations Number of iterations in each system.
	 * @param omega attenuation factor influencing convergence speed.
	 * @param norm The norm factor of the extended algebraic distance.
	 * @param withEdgeScores	 calculate array of scores for edges {u,v} that equal ad(u,v)
	 */
	AlgebraicDistance(const Graph& G, count numberSystems=10, count numberIterations=30, double omega=0.5, index norm=0, bool withEdgeScores=false);

	/**
	 *
	 */
	virtual void preprocess();

	/**
	 * @return algebraic distance between the two nodes.
	 */
	virtual double distance(node u, node v);


	virtual std::vector<double> getEdgeScores();


protected:

	/**
	 * initialize vectors randomly
	 */
	void randomInit();

	count numSystems; //!< number of vectors/systems used for algebraic iteration
	count numIters; //!< number of iterations in each system
	double omega; //!< attenuation factor influencing the speed of convergence
	index norm;
	const index MAX_NORM = 0;
	bool withEdgeScores;

	std::vector<double> loads; //!< loads[u*numSystems..(u+1)*numSystems]: loads for node u

	std::vector<double> edgeScores; //!< distance(u,v) for edge {u,v}

};

} /* namespace NetworKit */
#endif /* ALGEBRAICDISTANCE_H_ */
