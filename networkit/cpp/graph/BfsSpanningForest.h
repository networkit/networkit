/*
 * BfsSpanningForest.h
 *
 *  Created on: Aug 7, 2014
 *      Author: Christian Staudt
 */

#ifndef BFSSPANNINGFOREST_H_
#define BFSSPANNINGFOREST_H_

#include "Graph.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Creates a spanning forest (or tree).
 */
class BfsSpanningForest { // TODO: derive from SpanningForest when fixed
public:

	BfsSpanningForest(const Graph& G);

	Graph generate();

protected:

	const Graph& G;

};


} /* namespace NetworKit */
#endif /* BFSSPANNINGFOREST_H_ */
