/*
 * Sfigality.h
 *
 *  Created on: 20.01.2016
 *      Author:Christian Staudt
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
	 * Constructs the Sfigality class for the given Graph @a G. Sfigality indicates the number of neighbors with a higher degree than the node itself.
	 *
	 * @param G The graph.

	 */
	Sfigality(const Graph& G);

	void run() override;

	/**
	 * @return the theoretical maximum degree centrality, which is $n - 1$ (including the possibility of a self-loop)
	 */
	double maximum() override;
};

} /* namespace NetworKit */

#endif /* SFIGALITY_H_ */
