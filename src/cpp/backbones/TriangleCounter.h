/*
 * TriangleCounter.h
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner
 */

#ifndef TRIANGLECOUNTER_H_
#define TRIANGLECOUNTER_H_

#include "../graph/Graph.h"

namespace NetworKit {

/** 
 * Abstract base class for per-edge triangle counting algorithms.
 */
class TriangleCounter {

public:
	TriangleCounter();
	virtual ~TriangleCounter();

	/**
	 * Calculates triangle counts for the edges of the given graph.
	 */
	virtual std::vector<count> triangleCounts(const Graph& g) = 0;
};

} /* namespace NetworKit */
#endif /* TRIANGLECOUNTER_H_ */
