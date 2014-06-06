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

class EPPFactory {
public:
	EPP make(count ensembleSize, std::string baseAlgorithm, std::string finalAlgorithm);

};

}


#endif