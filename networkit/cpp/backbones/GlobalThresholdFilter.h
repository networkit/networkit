/*
 * GlobalThresholdFilter.h
 *
 *  Created on: 23.07.2014
 *      Author: Gerd Lindner
 */

#ifndef GLOBALTHRESHOLDFILTER_H_
#define GLOBALTHRESHOLDFILTER_H_

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Calculates a backbone by applying a global threshold to an edge attribute.
 */
class GlobalThresholdFilter {

public:

	/**
	 * Creates a new instance of a global threshold filter.
	 * @param threshold		the threshold
	 * @param above			if set to true, edge attribute needs to be above or equal to the threshold.
	 * 						If set to false, edge attribute needs to be below or equal to the threshold.
	 */
	GlobalThresholdFilter(double threshold, bool above); //TODO: better name for parameter?

	Graph calculate(const Graph& graph, const std::vector<double>& attribute);

private:
	double threshold;
	bool above;

	/**
	* Creates a new undirected graph that contains only the nodes of the given graph.
	* TODO: Implement a clone method in Graph instead?
	* @param graph 	the original graph to copy
	* @param weighted	whether the new graph should be weighted
	*/
	Graph cloneNodes(const Graph& graph, bool weighted);

};

}
/* namespace NetworKit */
#endif /* GLOBALTHRESHOLDFILTER_H_ */
