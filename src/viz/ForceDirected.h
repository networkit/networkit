/*
 * ForceDirected.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning
 */

#ifndef FORCEDIRECTED_H_
#define FORCEDIRECTED_H_

#include "SpringEmbedder.h"

namespace NetworKit {

class ForceDirected: public NetworKit::SpringEmbedder {
public:
	ForceDirected();
	virtual ~ForceDirected();

	/**
	 * Assigns coordinates to vertices in graph @a g
	 */
	virtual void draw(Graph& g);
};

} /* namespace NetworKit */
#endif /* FORCEDIRECTED_H_ */
