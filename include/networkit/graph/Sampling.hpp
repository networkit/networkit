/*
 * Sampling.h
 *
 *  Created on: 17.02.2014
 *      Author: cls
 */

#ifndef NETWORKIT_GRAPH_SAMPLING_HPP_
#define NETWORKIT_GRAPH_SAMPLING_HPP_

#include <networkit/graph/Graph.hpp>

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

#endif // NETWORKIT_GRAPH_SAMPLING_HPP_
