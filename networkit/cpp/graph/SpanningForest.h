/*
 * SpanningForest.h
 *
 *  Created on: 06.09.2015
 *      Author: Henning
 */

#ifndef SPANNINGFOREST_H_
#define SPANNINGFOREST_H_

#include "Graph.h"

namespace NetworKit {

/**
 * Base class for spanning forest/tree algorithms.
 */
class SpanningForest {
protected:
	const Graph& G;
	Graph forest;

public:
	SpanningForest(const Graph& G);
	virtual ~SpanningForest() = default;

	virtual void run();

	/**
	 * Deprecated. Please integrate into run method.
	 */
	Graph generate();

	/**
	 * @return Forest computed by run method.
	 * Note: So far no explicit check if run method has been invoked before.
	 */
	Graph getForest();
};

} /* namespace NetworKit */
#endif /* SPANNINGFOREST_H_ */
