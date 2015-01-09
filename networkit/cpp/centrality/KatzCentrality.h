/*
 * KatzCentrality.h
 *
 *  Created on: 09.01.2015
 *      Author: Henning
 */

#ifndef KATZCENTRALITY_H_
#define KATZCENTRALITY_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 * Computes the Katz centrality of the graph.
 */
class KatzCentrality: public Centrality {
protected:
	double alpha; // damping
	double beta; // constant centrality amount
	double tol; // error tolerance

public:
	/**
	 * Constructs the EigenvectorCentrality class for the given Graph @a G. @a tol defines the tolerance for convergence.
	 *
	 * @param[in] G The graph.
	 * @param[in] tol The tolerance for convergence.
	 */
	KatzCentrality(const Graph& G, double alpha = 1e-3, double beta = 1.0, double tol = 1e-9);

	virtual void run();
};

} /* namespace NetworKit */
#endif /* KATZCENTRALITY_H_ */
