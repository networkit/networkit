/*
 * SpanningEdgeCentrality.h
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
 * SpanningEdgeCentrality edge centrality.
 *
 */
class SpanningEdgeCentrality: public NetworKit::Centrality {
protected:
	double tol;
	Lamg lamg;
	uint64_t setupTime;

public:
	/**
	 * Constructs the SpanningEdgeCentrality class for the given Graph @a G.
	 * @param G The graph.
	 * @param tol constant used for the approximation: with probability at least 1-1/n, the approximated scores are within a factor 1+tol from the exact scores
	 */
	SpanningEdgeCentrality(const Graph& G, double tol = 0.1);

	/**
	 * Destructor.
	 */
	virtual ~SpanningEdgeCentrality() = default;


	/**
	* Compute spanning edge centrality scores exactly for all edges.
	*/
	void run() override;


	/**
	 * Compute approximation by JL projection.
	 */
	void runApproximation();

	/**
	 * Compute approximation by JL projection in parallel.
	 */
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
	/**
	 * Compute value for one edge only.
	 * @param[in] u Endpoint of edge.
	 * @param[in] v Endpoint of edge.
	 */
	double runForEdge(node u, node v);

};

} /* namespace NetworKit */


#endif /* SPANNING_H_ */
