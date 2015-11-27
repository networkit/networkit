/*
 * DynamicBarabasiAlbertGenerator.h
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef DYNAMICBARABASIALBERTGENERATOR_H_
#define DYNAMICBARABASIALBERTGENERATOR_H_

#include <set>

#include "DynamicGraphSource.h"

namespace NetworKit {


// FIXME: for k=2, degree 2 nodes should be most frequent but degree 4 nodes are
/**
 * @ingroup generators
 */
class DynamicBarabasiAlbertGenerator: public NetworKit::DynamicGraphSource {

protected:

	count k; 						//!< parameter of the BA model: number of edges per new node
	count degSum; 					//!< degree sum of current graph


public:


	DynamicBarabasiAlbertGenerator(count k = 2);

	virtual void initializeGraph();

	virtual void generate();

};

} /* namespace NetworKit */
#endif /* DYNAMICBARABASIALBERTGENERATOR_H_ */
