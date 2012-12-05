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

namespace EnsembleClustering {

class Generator {

public:

	Generator();

	virtual ~Generator();

	/**
	 * Generate a random graph according to the Erdšs-Renyi model.
	 *
	 * @param[in]	n	number of nodes
	 * @param[in]	p	edge probability
	 */
	Graph& makeErdosRenyiGraph(int64_t n, double p);


	/**
	 * Gemerate a graph whose nodes and edges form a circle.
	 *
	 * @param[in]	n	number of nodes
	 */
	Graph& makeCircularGraph(int64_t n);

};

} /* namespace EnsembleClustering */
#endif /* GENERATOR_H_ */
