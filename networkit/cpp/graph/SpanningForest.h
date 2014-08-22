/*
 * SpanningForest.h
 *
 *  Created on: Aug 7, 2014
 *      Author: Christian Staudt
 */

#ifndef SPANNINGFOREST_H_
#define SPANNINGFOREST_H_

#include "Graph.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Creates a spanning forest (or tree).
 */
class SpanningForest {
public:

	SpanningForest(const Graph& G);

	Graph generate();

	/** only to be used by cython - this eliminates an unnecessary copy */
	Graph* _generate() {
		return new Graph{std::move(generate())};
	};

protected:

	const Graph& G;

};


} /* namespace NetworKit */
#endif /* SPANNINGFOREST_H_ */
