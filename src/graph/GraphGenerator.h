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
#include "../clustering/Clustering.h"

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
	Graph makeErdosRenyiGraph(int64_t n, double p);


	/**
	 * Generate a graph whose nodes and edges form a circle.
	 *
	 * @param[in]	n	number of nodes
	 */
	Graph makeCircularGraph(int64_t n);


	/**
	 * Generate a complete graph.
	 *
	 * @param[in]	n	number of nodes
	 */
	Graph makeCompleteGraph(int64_t n);


	/**
	 * Creates a clustered random graph:
	 *
	 * @param[in]	n	number of nodes
	 * @param[in]	k	number of clusters
	 * @param[in]	pin		intra-cluster edge probability
	 * @param[in]	pout	inter-cluster edge probability
	 */
	Graph makeClusteredRandomGraph(int64_t n, int64_t k, double pin, double pout);

};

} /* namespace EnsembleClustering */
#endif /* GENERATOR_H_ */
