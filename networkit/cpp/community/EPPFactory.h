/*
 * EPPFactory.h
 *
 *  Created on: 12.11.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef EPPFACTORY_H_
#define EPPFACTORY_H_

#include "EPP.h"
namespace NetworKit {
namespace EPPFactory {
	/**
	 * Creates an instance of the ensemble preprocessing.
	 * @param	G	The graph on which the ensemble is supposed to run.
	 * @param	ensembleSize	The amount of baseAlgorithms to preprocess the communities.
	 * @param	baseAlgorithm	String representation of the algorithm ("PLP","PLM") to preprocess the communities. ensembleSize instances will be created.
	 * @param	finalAlgorithm	String representation of the algorithm ("PLP" "PLM[R]") to finish the ensemble.
	 */
	EPP make(const Graph& G, count ensembleSize, std::string baseAlgorithm, std::string finalAlgorithm);

	/**
	 * Creates an instance of the ensemble preprocessing and returns a pointer to it (for Cython).
	 * @param	G	The graph on which the ensemble is supposed to run.
	 * @param	ensembleSize	The amount of baseAlgorithms to preprocess the communities.
	 * @param	baseAlgorithm	String representation of the algorithm ("PLP","PLM") to preprocess the communities. ensembleSize instances will be created.
	 * @param	finalAlgorithm	String representation of the algorithm ("PLP" "PLM[R]") to finish the ensemble.
	 */
	EPP* makePtr(const Graph& G, count ensembleSize, std::string baseAlgorithm, std::string finalAlgorithm);
}
}
#endif 