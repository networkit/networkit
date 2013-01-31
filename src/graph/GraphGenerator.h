/*
 * Generator.h
 *
 *  Created on: 05.12.2012
 *      Author: cls
 */

#ifndef GENERATOR_H_
#define GENERATOR_H_

#include "Graph.h"
#include "../aux/RandomProbability.h"
#include "../aux/RandomInteger.h"
#include "../clustering/base/Clustering.h"

namespace EnsembleClustering {

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
	virtual Graph makeErdosRenyiGraph(int64_t n, double p);


	/**
	 * Alias for makeErdosRenyiGraph
	 */
	virtual Graph makeRandomGraph(int64_t n, double p);


	/**
	 * Generate a graph whose nodes and edges form a circle.
	 *
	 * @param[in]	n	number of nodes
	 */
	virtual Graph makeCircularGraph(int64_t n);


	/**
	 * Generate a complete graph.
	 *
	 * @param[in]	n	number of nodes
	 */
	virtual Graph makeCompleteGraph(int64_t n);


	/**
	 * Creates a clustered random graph:
	 *
	 * @param[in]	n	number of nodes
	 * @param[in]	k	number of clusters
	 * @param[in]	pin		intra-cluster edge probability
	 * @param[in]	pout	inter-cluster edge probability
	 */
	virtual Graph makeClusteredRandomGraph(int64_t n, int64_t k, double pin, double pout);

	/**
	 * Creates a clustered random graph:
	 *
	 * @param[in]	n	number of nodes
	 * @param[in]	k	number of clusters
	 * @param[in]	pin		intra-cluster edge probability
	 * @param[in]	pout	inter-cluster edge probability
	 */
	virtual std::pair<Graph, Clustering> makeClusteredRandomGraphWithReferenceClustering(int64_t n, int64_t k, double pin, double pout);


	/**
	 * Create a clustered random graph from a given clustering.
	 *
	 */
	virtual Graph makeClusteredRandomGraph(Clustering& zeta, double pin, double pout);



	/**
	 * Generate random graph according to the Barabasi-Albert model (preferential attachment)
	 *
	 *
	 */
	virtual Graph makeBarabasiAlbertGraph(int64_t n, int64_t k);

};

} /* namespace EnsembleClustering */
#endif /* GENERATOR_H_ */
