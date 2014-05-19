/*
 * Sampling.h
 *
 *  Created on: 17.02.2014
 *      Author: cls
 */

#ifndef SAMPLING_H_
#define SAMPLING_H_

 #include "Graph.h"

namespace NetworKit {

class Sampling {

public:

	static node randomNode(const IGraph& G);

	static std::pair<node, node> randomEdge(const IGraph& G);

	static node randomNeighbor(const IGraph& G, node u);

};

} /* namespace NetworKit */

#endif /* SAMPLING_H_ */
