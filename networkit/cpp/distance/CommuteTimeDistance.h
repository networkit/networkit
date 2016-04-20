/*
 * CommuteTimeDistance.h
 *
 *  Created on: 12.04.2016
 *      Author: ebergamini
 */

#ifndef COMMUTETIMEDIST_H_
#define COMMUTETIMEDIST_H_

#include "../numerics/LAMG/Lamg.h"
#include "../graph/Graph.h"
#include "../base/Algorithm.h"


namespace NetworKit {

/**
 * @ingroup centrality
 *
 * CommuteTimeDistance edge centrality.
 *
 */
class CommuteTimeDistance: public Algorithm {

public:
	/**
	 * Constructs the CommuteTimeDistance class for the given Graph @a G.
	 * @param G The graph.
	 */
	CommuteTimeDistance(const Graph& G, double tol = 0.1);

	/**
	 * Destructor.
	 */
	virtual ~CommuteTimeDistance() = default;


	/**
	 * Computes ECTD exactly.
	 */
	virtual void run();
	/**
	 * Compute approximation by projection.
	 */
	void runApproximation();

	void runParallelApproximation();
	/**
	 * @return The elapsed time to setup the solver in milliseconds.
	 */
	uint64_t getSetupTime() const;
	/**
	 * Returns the commute time distance between node @a u and node @a v.
	 * @return commute time distance between the two nodes.
	 */
 	double distance(node u, node v);
	/**
	 * Returns the commute time distance between node @a u and node @a v.
	 This method does not need the initial preprocessing.
	 * @return commute time distance between the two nodes.
	 */
	double runSinglePair(node u, node v);

protected:
	const Graph& G;
	double tol;
	Lamg lamg;
	uint64_t setupTime;
	std::vector<std::vector<double>> distances;
	std::vector<Vector> solutions;
	bool hasRun = false;
	bool exactly;
	count k;
};

} /* namespace NetworKit */


#endif /* COMMUTETIMEDIST_H_ */
