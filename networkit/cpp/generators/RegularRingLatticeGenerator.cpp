/*
* RegularRingLatticeGenerator.cpp
*
*  Created on: 09.07.2014
*      Author: Simon Bischof
*/

#include "RegularRingLatticeGenerator.h"

namespace NetworKit {

RegularRingLatticeGenerator::RegularRingLatticeGenerator(count nNodes, count nNeighbors)
	: nNodes(nNodes), nNeighbors(nNeighbors) {
	if (nNeighbors >= nNodes / 2 - 1) {
		nNeighbors = nNodes / 2 - 1;
	}
}


Graph RegularRingLatticeGenerator::generate() {
	Graph G(nNodes);
	for (count i = 0; i < nNodes; i++) {
		for (count j = 1; j <= nNeighbors; j++) {
			G.addEdge(i, (i + j) % nNodes);
		}
	}
	
	return G;
}

} /* namespace NetworKit */
