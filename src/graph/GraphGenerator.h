/*
 * Generator.h
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GENERATOR_H_
#define GENERATOR_H_

#include "Graph.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 * Simple methods for graph generation.
 */
class GraphGenerator {

public:

	GraphGenerator();

	virtual ~GraphGenerator();

	/**
	 * Generate a random graph according to the Erdos-Renyi model.
	 *
	 * @param[in]	n	number of nodes
	 * @param[in]	p	edge probability
	 */
	virtual Graph makeErdosRenyiGraph(count n, double p);


	/**
	 * Alias for makeErdosRenyiGraph
	 */
	virtual Graph makeRandomGraph(count n, double p);


	/**
	 * Generate a graph whose nodes and edges form a circle.
	 *
	 * @param[in]	n	number of nodes
	 */
	virtual Graph makeCircularGraph(count n);


	/**
	 * Generate a complete graph.
	 *
	 * @param[in]	n	number of nodes
	 */
	virtual Graph makeCompleteGraph(count n);


	/**
	 * Creates a clustered random graph:
	 *
	 * @param[in]	n	number of nodes
	 * @param[in]	k	number of clusters
	 * @param[in]	pin		intra-cluster edge probability
	 * @param[in]	pout	inter-cluster edge probability
	 */
	virtual Graph makeClusteredRandomGraph(count n, count k, double pin, double pout);

	/**
	 * Creates a clustered random graph:
	 *
	 * @param[in]	n	number of nodes
	 * @param[in]	k	number of clusters
	 * @param[in]	pin		intra-cluster edge probability
	 * @param[in]	pout	inter-cluster edge probability
	 */
	virtual std::pair<Graph, Partition> makeClusteredRandomGraphWithReferenceClustering(count n, count k, double pin, double pout);


	/**
	 * Create a clustered random graph from a given clustering.
	 *
	 */
	virtual Graph makeClusteredRandomGraph(Partition& zeta, double pin, double pout);


};

} /* namespace NetworKit */
#endif /* GENERATOR_H_ */
