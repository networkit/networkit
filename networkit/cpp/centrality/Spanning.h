/*
 * Spanning.h
 *
 *  Created on: 29.07.2015
 *      Author: henningm
 */

#ifndef SPANNING_H_
#define SPANNING_H_

#include "Centrality.h"
#include "../numerics/LAMG/SolverLamg.h"


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
	Smoother* smoother;
	SolverLamg* solver;
	LevelHierarchy* hierarchy;

public:
	/**
	 * Constructs the Spanning class for the given Graph @a G.
	 * @param G The graph.
	 */
	Spanning(const Graph& G, double tol = 1e-8);

	/**
	 * Destructor.
	 */
	virtual ~Spanning();


	/**
	* Compute spanning scores for all edges to desired tolerance.
	*/
	void run() override;

	/**
	 * Compute approximation by projection.
	 */
	void runApproximation();

	void runTreeApproximation();

	void runPseudoTreeApproximation();

	/**
	 * Compute value for one edge only.
	 * @param[in] u Endpoint of edge.
	 * @param[in] v Endpoint of edge.
	 */
	double runForEdge(node u, node v);

};

} /* namespace NetworKit */


#endif /* SPANNING_H_ */
