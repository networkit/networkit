/*
 * GraphProperties.h
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#ifndef GRAPHPROPERTIES_H_
#define GRAPHPROPERTIES_H_


#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup properties
 * @deprecated DEPRECATED: Implement algorithms in their own classes.
 * 
 * Collection of methods for basic network properties.
 */
class GraphProperties {
protected:
	/**
	 * @return Degree assortativity of the graph @a G.
	 * Based on Eq. (4) in Newman: Assortative mixing in networks.
	 * URL: http://arxiv.org/pdf/cond-mat/0205405.pdf.
	 * A similar description of this algorithm can be found in
	 * Newman: Networks. An Introduction. Chapter 8.7.
	 */
	static double degreeAssortativitySlower(const Graph& G, bool useWeights = false);


public:
	GraphProperties();
	virtual ~GraphProperties();

	static std::vector<count> degreeDistribution(const Graph& G);

	static std::vector<count> degreeSequence(const Graph& G);


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
	static std::vector<double> localClusteringCoefficients(const Graph& G);


	/**
	 * The average local clustering coefficient for the graph @a G.
	 * 		$\frac{1}{n} \cdot \sum_{v \in V} c_v$
	 *
	 * @param[in]	G	the graph
	 * @return Average local clustering coefficient.
	 */
	static double averageLocalClusteringCoefficient(const Graph& G);

	static std::vector<double> localClusteringCoefficientPerDegree(const Graph& G);

	static std::pair<count, count> minMaxDegree(const Graph& G);

	static std::pair<std::pair<count,count>, std::pair<count,count>> minMaxDegreeDirected(const Graph& G);

	static double averageDegree(const Graph& G);

	/**
	 * Get degree assortativity of the graph @a G.
	 *
	 * @param G The graph
	 * @param useWeights If @c true, the weights are considered for calculation. Default: @c false.
	 * @return Degree assortativity of the graph @a G.
	 * @note Degree assortativity based on description in Newman: Networks. An Introduction. Chapter 8.7.
	 */
	static double degreeAssortativity(const Graph& G, bool useWeights = false);


	/**
	 * Degree Assortativity for directed graphs.[1]
	 *
	 * [1] http://www.pnas.org/content/107/24/10815.full
	 * @param  G     The graph
	 * @param  alpha True, if the out-degree of the edge's source vertex is to be considered. False for the in-degree.
	 * @param  beta  True, if the out-degree of the edge's target vertex is to be considered. False for the in-degree.
	 * @return       The degree assortativity for the given configuration.
	 */
	static double degreeAssortativityDirected(const Graph& G, bool alpha, bool beta);
};

} /* namespace NetworKit */
#endif /* GRAPHPROPERTIES_H_ */
