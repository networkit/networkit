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

/**
 * @ingroup graph
 */
class Sampling {

public:

	static node randomNode(const Graph& G);

	static std::pair<node, node> randomEdge(const Graph& G);

	static node randomNeighbor(const Graph& G, node u);

};

} /* namespace NetworKit */

#endif /* SAMPLING_H_ */
