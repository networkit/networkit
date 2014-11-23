/*
 * EnsembleFactory.cpp
 *
 *  Created on: 12.11.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef EPPFACTORY_H_
#define EPPFACTORY_H_

#include "EPP.h"
namespace NetworKit {
namespace EPPFactory {
	EPP make(const Graph& G, count ensembleSize, std::string baseAlgorithm, std::string finalAlgorithm);
	EPP* makePtr(const Graph& G, count ensembleSize, std::string baseAlgorithm, std::string finalAlgorithm);
}
}
#endif 