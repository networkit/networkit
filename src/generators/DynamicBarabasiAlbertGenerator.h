/*
 * DynamicBarabasiAlbertGenerator.h
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef DYNAMICBARABASIALBERTGENERATOR_H_
#define DYNAMICBARABASIALBERTGENERATOR_H_

#include <set>

#include "DynamicGraphGenerator.h"
#include "../aux/RandomInteger.h"

namespace NetworKit {

class DynamicBarabasiAlbertGenerator: public NetworKit::DynamicGraphGenerator {

protected:

	Aux::RandomInteger randInt;		//!< random integer generators
	count k; 						//!< parameter of the BA model: number of edges per new node
	count degSum; 					//!< degree sum of current graph


public:

	DynamicBarabasiAlbertGenerator(GraphEventProxy& proxy, count k = 2);

	virtual ~DynamicBarabasiAlbertGenerator();

	virtual void initializeGraph();

	virtual void generateWhile(std::function<bool(void)> cont);
};

} /* namespace NetworKit */
#endif /* DYNAMICBARABASIALBERTGENERATOR_H_ */
