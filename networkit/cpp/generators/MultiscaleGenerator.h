/*
 * MultiscaleGenerator.h
 *
 *  Created on: Oct 27, 2015
 *      Author: Christian Staut
 */

#ifndef MULTISCALEGENERATOR_H_
#define MULTISCALEGENERATOR_H_


#include "StaticGraphGenerator.h"

namespace NetworKit {

/**
 * @ingroup generators
 * TODO: doc
 */
class MultiscaleGenerator: public NetworKit::StaticGraphGenerator {
private:

public:
	MultiscaleGenerator();

	Graph generate() override;
};

} /* namespace NetworKit */
#endif /* MULTISCALEGENERATOR */
