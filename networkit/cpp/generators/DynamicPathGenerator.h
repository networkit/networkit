/*
 * DynamicPathGenerator.h
 *
 *  Created on: 14.01.2014
 *      Author: cls
 */

#ifndef DYNAMICPATHGENERATOR_H_
#define DYNAMICPATHGENERATOR_H_

#include "DynamicGraphGenerator.h"
#include <map>

namespace NetworKit {

/**
 * @ingroup generators
 * Example dynamic graph generator: Generates a dynamically growing path.
 */
class DynamicPathGenerator: public NetworKit::DynamicGraphGenerator {
public:

	std::vector<GraphEvent> generate(count nSteps) override;

};

} /* namespace NetworKit */

#endif /* DYNAMICPATHGENERATOR_H_ */
