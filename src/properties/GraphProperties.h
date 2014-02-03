/*
 * GraphProperties.h
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#ifndef GRAPHPROPERTIES_H_
#define GRAPHPROPERTIES_H_


#include "../graph/Graph.h"
#include "../io/METISGraphReader.h"

namespace NetworKit {

/**
 * Collection of methods for basic network properties.
 */
class GraphProperties {
public:
	GraphProperties();
	virtual ~GraphProperties();

	static std::vector<count> degreeDistribution(Graph& G);

	static std::vector<unsigned int> degreeSequence(Graph& G); // TODO: revert to count when cython issue fixed


	/**
	 * The local clustering coefficient for a node is the number of edges among its
	 * neighbors divided by the number of possible edges.
	 * 	For an undirected graph where N(v) does not include v itself:
	 * 		$c_v := \frac{2 \cdot |E(N(v))| }{\deg(v) ( \deg(v) - 1 )}$
	 *
	 *
	 * @param[in]	G	the graph
	 * @param[out]		node -> local clustering coefficient
	 */
	static std::vector<double> localClusteringCoefficients(Graph& G);


	/**
	 * The average local clustering coefficient for the graph.
	 * 		$\frac{1}{n} \cdot \sum_{v \in V} c_v$
	 *
	 * @param[in]	G	the graph
	 */
	static double averageLocalClusteringCoefficient(Graph& G);

	static std::vector<double> localClusteringCoefficientPerDegree(Graph& G);

	static std::pair<count, count> minMaxDegree(Graph& G);

	static double averageDegree(const Graph& G);
};

} /* namespace NetworKit */
#endif /* GRAPHPROPERTIES_H_ */
