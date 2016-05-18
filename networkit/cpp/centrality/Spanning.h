/*
 * Spanning.h
 *
 *  Created on: 29.07.2015
 *      Author: henningm
 */

#ifndef SPANNING_H_
#define SPANNING_H_

#include "Centrality.h"
#include "../numerics/LAMG/Lamg.h"


namespace NetworKit {

/**
 * @ingroup centrality
 *
 * Spanning edge centrality.
 *
 */
class Spanning: public NetworKit::Centrality {
protected:
	double tol;
	Lamg lamg;
	uint64_t setupTime;

public:
	/**
	 * Constructs the Spanning class for the given Graph @a G.
	 * @param G The graph.
	 */
	Spanning(const Graph& G, double tol = 0.1);

	/**
	 * Destructor.
	 */
	virtual ~Spanning() = default;


	/**
	* Compute spanning scores for all edges to desired tolerance.
	*/
	void run() override;


	/**
	 * Compute approximation by projection.
	 */
	void runApproximation();

	void runParallelApproximation();

	/**
	 * Only used by benchmarking. Computes an approximation by projection and solving Laplacian systems.
	 * Measures the time needed to compute the approximation and writes the problem vectors to the
	 * directory of the graph specified by @a graphPath.
	 * @param directory
	 * @return Elapsed time in milliseconds.
	 */
	uint64_t runApproximationAndWriteVectors(const std::string &graphPath);

	/**
	 * @return The elapsed time to setup the solver in milliseconds.
	 */
	uint64_t getSetupTime() const;

	void runTreeApproximation(count reps);

	void runTreeApproximation2(count reps);

	void runPseudoTreeApproximation(count reps);

	/**
	 * Compute value for one edge only.
	 * @param[in] u Endpoint of edge.
	 * @param[in] v Endpoint of edge.
	 */
	double runForEdge(node u, node v);

};

} /* namespace NetworKit */


#endif /* SPANNING_H_ */
