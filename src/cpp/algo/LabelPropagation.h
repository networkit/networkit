/*
* LabelPropagation.h
*
* A generic implementation of parallel label propagation.
*
*  Created on: 04.08.2014
*      Author: Christian Staudt (christian.staudt@kit.edu)
*/

#ifndef LABELPROPAGATION_H_
#define LABELPROPAGATION_H_

#include "../structures/Partition.h"


namespace NetworKit {

typedef index label;

class LabelPropagation {

protected:

	std::vector<label> labeling;	//!< a labelling of the node set
	Graph G;

public:

	void run();
};

}

#endif /* LABELPROPAGATION_H_ */
