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
 * Calculates a sparsified graph by applying a global threshold to an edge score.
 */
class GlobalThresholdFilter {

public:

	/**
	 * Creates a new instance of a global threshold filter.
	 * @param threshold		the threshold
	 * @param above			if set to true, edges with a score above or equal to the threshold remain in the filtered graph.
	 * 						If set to false, edges with a score below or equal to the threshold remain in the filtered graph.
	 */
	GlobalThresholdFilter(const Graph& graph, const std::vector<double>& attribute, double threshold, bool above);

	Graph calculate();

private:
	const Graph& graph;
	const std::vector<double>& attribute;
	double threshold;
	bool above;

	/**
	* Creates a new undirected graph that contains only the nodes of the given graph.
	* @param graph 	the original graph to copy
	* @param weighted	whether the new graph should be weighted
	*/
	Graph cloneNodes(const Graph& graph, bool weighted);

};

}
/* namespace NetworKit */
#endif /* GLOBALTHRESHOLDFILTER_H_ */
