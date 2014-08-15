/*
* WattsStrogatzGenerator.h
*
*  Created on: 09.07.2014
*      Author: Simon Bischof
*/

#ifndef WattsStrogatzGENERATOR_H_
#define WattsStrogatzGENERATOR_H_

#include "StaticGraphGenerator.h"

namespace NetworKit {


class WattsStrogatzGenerator: public StaticGraphGenerator {

public:
	/**
	* Constructs a graph according to the Watts and Strogatz model
	* (https://en.wikipedia.org/wiki/Watts_and_Strogatz_model),
	* which produces graphs with high clustering and low average path length.
	*
	* First, a regular ring lattice is generated.
	* Then some edges are rewired randomly.
	*
	* @param nNodes 		number of nodes in target graph
	* @param nNeighbors		number of neighbors on each side of a node
	* @param p				rewiring probability
	*/
	WattsStrogatzGenerator(count nNodes, count nNeighbors, double p);

	virtual Graph generate();

protected:
		count nNodes;
		count nNeighbors;
		double p;
};

} /* namespace NetworKit */
#endif /* WattsStrogatzGENERATOR_H_ */
