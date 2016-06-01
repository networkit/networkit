/*
 * Sfigality.h
 *
 *  Created on: 20.01.2016
 *      Author: Elisabetta Bergamini, Christian Staudt
 */

#ifndef SFIGALITY_H_
#define SFIGALITY_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 * A
 */
class Sfigality: public NetworKit::Centrality {
public:
	/**
	 * Constructs the Sfigality class for the given Graph @a G.
	 *
	 * @param G The graph.

	 */
	Sfigality(const Graph& G);

	void run() override;

	/**
	 * @return the theoretical maximum degree centrality, which is $n$ (including the possibility of a self-loop)
	 */
	double maximum() override;
};

} /* namespace NetworKit */

#endif /* SFIGALITY_H_ */
