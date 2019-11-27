/*
 * Sampling.cpp
 *
 *  Created on: 17.02.2014
 *      Author: cls
 */

#include <networkit/graph/Sampling.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/auxiliary/Random.hpp>

namespace NetworKit {

node Sampling::randomNode(const Graph& G) {
    return GraphTools::randomNode(G);
}

// the following methdods are commented in order to create linker-errors should they be used before
// they are actually defined (not returning from a function with a returntype != void is UB):

//std::pair<node, node> Sampling::randomEdge(const Graph& G) {
//	// TODO: implement
//}

//node Sampling::randomNeighbor(const Graph& G, node u) {
//	// TODO: implement
//}



} /* namespace NetworKit */

