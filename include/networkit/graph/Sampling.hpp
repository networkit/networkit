/*
 * Sampling.hpp
 *
 *  Created on: 17.02.2014
 *      Author: cls
 */

#ifndef NETWORKIT_GRAPH_SAMPLING_HPP_
#define NETWORKIT_GRAPH_SAMPLING_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

/**
 * @ingroup graph
 */
class Sampling {

public:

    static node randomNode(const Graph& G) {
        return GraphTools::randomNode(G);
    }

    static std::pair<node, node> randomEdge(const Graph& G);
        return GraphTools::randomEdge(G);
    }

    static node randomNeighbor(const Graph& G, node u);
        return GraphTools::randomNeighbor(G, u);
    }

};

} /* namespace NetworKit */

#endif // NETWORKIT_GRAPH_SAMPLING_HPP_
