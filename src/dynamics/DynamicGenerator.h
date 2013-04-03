/*
 * DynamicGenerator.h
 *
 *  Created on: 03.04.2013
 *      Author: cls
 */

#ifndef DYNAMICGENERATOR_H_
#define DYNAMICGENERATOR_H_

#include <functional>

namespace NetworKit {

class DynamicGenerator {

public:

	DynamicGenerator();

	virtual ~DynamicGenerator();

	virtual void generate(std::function<bool(void)> terminate) = 0;
};

} /* namespace NetworKit */
#endif /* DYNAMICGENERATOR_H_ */
