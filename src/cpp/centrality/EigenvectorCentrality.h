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
 * Computes the leading eigenvector of the graph's adjacency matrix (normalized in 2-norm).
 * Interpreted as eigenvector centrality score.
 */
class EigenvectorCentrality: public Centrality {
protected:
	double tol; // error tolerance

public:
	/**
	 * @param[in] G
	 * @param[in] normalized True if scores should be normalized in the interval [0,1].
	 */
	EigenvectorCentrality(const Graph& G, double tol = 1e-9);

	virtual void run();
};

} /* namespace NetworKit */
#endif /* EIGENVECTORCENTRALITY_H_ */
