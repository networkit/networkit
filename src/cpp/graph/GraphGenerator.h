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
 * @ingroup graph
 * Simple methods for graph generation.
 * DEPRECATED: do not use in new code, replace usage with dedicated generator classes
 */
class GraphGenerator {

public:

	/** Default destructor */
	virtual ~GraphGenerator() = default;

	/**
	 * Generate a random graph according to the Erdos-Renyi model.
	 *
	 * @param[in]	n	Number of nodes.
	 * @param[in]	p	Edge probability.
	 * @return The generated graph.
	 */
	virtual Graph makeErdosRenyiGraph(count n, double p);


	/**
	 * Alias for @link makeErdosRenyiGraph(count n, double p)
	 */
	virtual Graph makeRandomGraph(count n, double p);


	/**
	 * Generate a graph whose nodes and edges form a circle.
	 *
	 * @param[in]	n	Number of nodes.
	 * @return The generated graph.
	 */
	virtual Graph makeCircularGraph(count n);


	/**
	 * Generate a complete graph.
	 *
	 * @param[in]	n	Number of nodes.
	 * @return The generated graph.
	 */
	virtual Graph makeCompleteGraph(count n);


	/**
	 * Creates a clustered random graph.
	 *
	 * @param[in]	n	Number of nodes.
	 * @param[in]	k	Number of clusters.
	 * @param[in]	pin		Intra-cluster edge probability.
	 * @param[in]	pout	Inter-cluster edge probability.
	 * @return The generated graph.
	 */
	virtual Graph makeClusteredRandomGraph(count n, count k, double pin, double pout);

	/**
	 * Creates a clustered random graph and a reference clustering for this graph.
	 *
	 * @param[in]	n	Number of nodes.
	 * @param[in]	k	Number of clusters.
	 * @param[in]	pin		Intra-cluster edge probability.
	 * @param[in]	pout	Inter-cluster edge probability.
	 * @return The generated graph and a reference clustering of it.
	 */
	virtual std::pair<Graph, Partition> makeClusteredRandomGraphWithReferenceClustering(count n, count k, double pin, double pout);


	/**
	 * Create a clustered random graph from a given clustering @a zeta.
	 *
	 * @param[in] 	zeta	The clustering.
	 * @param[in] 	pin		Intra-cluster edge probability.
	 * @param[in]	pout	Inter-cluster edge probability.
	 * @return The generated graph.
	 */
	virtual Graph makeClusteredRandomGraph(Partition& zeta, double pin, double pout);


};

} /* namespace NetworKit */
#endif /* GENERATOR_H_ */
