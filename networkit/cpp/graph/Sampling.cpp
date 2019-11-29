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

} /* namespace NetworKit */

