/*
 * RandomSpanningForest.h
 *
 *  Created on: 06.09.2015
 *      Author: Henning
 */

#ifndef RANDOMSPANNINGFOREST_H_
#define RANDOMSPANNINGFOREST_H_

#include "../graph/Graph.h"
#include "../graph/SpanningForest.h"

namespace NetworKit {

/**
 * Creates a random spanning tree for each connected component.
 * Time complexity: cover time of G.
 * @ingroup graph
 */
class RandomSpanningForest: public SpanningForest {
public:
	RandomSpanningForest(const Graph& G);
	virtual ~RandomSpanningForest() = default;

	/**
	 * Computes for each component a random spanning tree.
	 * Uses simple random-walk based algorithm.
	 * Time complexity: cover time of G.
	 */
	virtual void run() override;
};

} /* namespace NetworKit */
#endif /* RANDOMSPANNINGFOREST_H_ */
