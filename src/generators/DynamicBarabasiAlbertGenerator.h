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
#include "../auxiliary/RandomInteger.h"

namespace NetworKit {


// FIXME: for k=2, degree 2 nodes should be most frequent but degree 4 nodes are
class DynamicBarabasiAlbertGenerator: public NetworKit::DynamicGraphGenerator {

protected:

	Aux::RandomInteger randInt;		//!< random integer generators
	count k; 						//!< parameter of the BA model: number of edges per new node
	count degSum; 					//!< degree sum of current graph


public:

	DynamicBarabasiAlbertGenerator(GraphEventProxy& proxy, count k = 2);

	virtual ~DynamicBarabasiAlbertGenerator();

	virtual void initializeGraph();

	virtual void generate();

};

} /* namespace NetworKit */
#endif /* DYNAMICBARABASIALBERTGENERATOR_H_ */
