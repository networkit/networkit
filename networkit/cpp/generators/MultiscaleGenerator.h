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

protected:

    // parameters
    const Graph& original;
    // const double nodeGrowthRate;
    // const double edgeGrowthRate;
    // const double nodeEditingRate;
    // const double edgeEditingRate;
    const count maxLevels = 2;     // TODO: levels of hierarchy
    std::string aggregationScheme = "matching";


public:
    // TODO: parameters: retainIntermediates, coarseningScheme,
	MultiscaleGenerator(const Graph& original);

	Graph generate() override;
};

} /* namespace NetworKit */
#endif /* MULTISCALEGENERATOR */
