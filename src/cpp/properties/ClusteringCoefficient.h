/*
 * ClusteringCoefficient.h
 *
 *  Created on: 08.04.2013
 *      Author: Lukas Barth, David Weiß
 */

#ifndef CLUSTERINGCOEFFICIENT_H_
#define CLUSTERINGCOEFFICIENT_H_

#include "../graph/Graph.h"

namespace NetworKit {

class ClusteringCoefficient {

public:

	static std::vector<double> exactLocal(Graph &G);
	
	/**
	 * This calculates the average local clustering coefficient
	 * $$c(G) := \frac{1}{n} \sum_{u \in V} c(u)$$
	 *
	 * where $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$
	 */
	static double avgLocal(Graph& G);
  	static double approxAvgLocal(Graph& G, const count trials);

	/**
	 * This calculates the global clustering coefficient
	 */
  	static double exactGlobal(Graph& G);
  	static double approxGlobal(Graph& G, const count trials);
  
};

} /* namespace NetworKit */
#endif /* CLUSTERINGCOEFFICIENT_H_ */
