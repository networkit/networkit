/*
 * EigenvectorCentrality.h
 *
 *  Created on: 19.03.2014
 *      Author: Henning
 */

#ifndef EIGENVECTORCENTRALITY_H_
#define EIGENVECTORCENTRALITY_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 * Computes the leading eigenvector of the graph's adjacency matrix (normalized in 2-norm).
 * Interpreted as eigenvector centrality score.
 */
class EigenvectorCentrality: public Centrality {
protected:
	double tol; // error tolerance

public:
	/**
	 * Constructs an EigenvectorCentrality object for the given Graph @a G. @a tol defines the tolerance for convergence.
	 *
	 * @param[in] G The graph.
	 * @param[in] tol The tolerance for convergence.
	 */
	EigenvectorCentrality(const Graph& G, double tol = 1e-8);

	virtual void run();
};

} /* namespace NetworKit */
#endif /* EIGENVECTORCENTRALITY_H_ */
