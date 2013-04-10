/*
 * DynamicBarabasiAlbertGenerator.h
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef DYNAMICBARABASIALBERTGENERATOR_H_
#define DYNAMICBARABASIALBERTGENERATOR_H_

#include "DynamicGraphGenerator.h"
#include "../aux/RandomInteger.h"

namespace NetworKit {

class DynamicBarabasiAlbertGenerator: public NetworKit::DynamicGraphGenerator {

protected:

	int k; //<! paramter of the BA model: number of edges per new node


public:

	DynamicBarabasiAlbertGenerator(GraphEventProxy& proxy, int k = 2);

	virtual ~DynamicBarabasiAlbertGenerator();

	virtual void initializeGraph();

	virtual void generate(std::function<bool(void)> terminate);
};

} /* namespace NetworKit */
#endif /* DYNAMICBARABASIALBERTGENERATOR_H_ */
