/*
 * DynamicBarabasiAlbertGenerator.h
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef DYNAMICBARABASIALBERTGENERATOR_H_
#define DYNAMICBARABASIALBERTGENERATOR_H_

#include "DynamicGenerator.h"

namespace NetworKit {

class DynamicBarabasiAlbertGenerator: public NetworKit::DynamicGraphGenerator {

protected:

	virtual void initializeGraph();

public:

	DynamicBarabasiAlbertGenerator();

	virtual ~DynamicBarabasiAlbertGenerator();

	virtual void generate(std::function<bool(void)> terminate);
};

} /* namespace NetworKit */
#endif /* DYNAMICBARABASIALBERTGENERATOR_H_ */
