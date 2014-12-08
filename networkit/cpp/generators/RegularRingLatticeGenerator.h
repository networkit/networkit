/*
* RegularRingLatticeGenerator.h
*
*  Created on: 09.07.2014
*      Author: Simon Bischof
*/

#ifndef REGULARRINGLATTICEGENERATOR_H_
#define REGULARRINGLATTICEGENERATOR_H_

#include "StaticGraphGenerator.h"

namespace NetworKit {


class RegularRingLatticeGenerator: public StaticGraphGenerator {

public:
	/**
	* Construct a undirected regular ring lattice.
	*
	* @param nNodes 		number of nodes in target graph
	* @param nNeighbors		number of neighbors on each side of a node
	*/
	RegularRingLatticeGenerator(count nNodes, count nNeighbors);

	virtual Graph generate();

protected:
		count nNodes;
		count nNeighbors;

};

} /* namespace NetworKit */
#endif /* REGULARRINGLATTICEGENERATOR_H_ */
