/*
 * Graph.hpp
 *
 *  Created on: 24.05.2024
 *      Author: Fabian Brandt-Tumescheit
 */

#ifndef NETWORKIT_GRAPH_HYPERGRAPH_HPP_
#define NETWORKIT_GRAPH_HYPERGRAPH_HPP_

#include <networkit/Globals.hpp>
#include <networkit/graph/Attributes.hpp>
#include <networkit/graph/EdgeIterators.hpp>
#include <networkit/graph/NeighborIterators.hpp>
#include <networkit/graph/NodeIterators.hpp>

namespace NetworKit {

/**
 * @ingroup graph
 * A hypergraph (with optional weights) and parallel iterator methods.
 */
class Hypergraph final {

    // Hypergraph basic attributes

    // current number of nodes
    count n;

    // current number of edges
    count m;

    //!< current upper bound of node ids, z will be the id of the next node
    node z;
    //!< current upper bound of edge ids, will be the id of the next edge
    edgeid omega;

    // Node related data

    //!< nodeExists[v] is false if node v has been removed from the graph
    std::vector<bool> nodeExists;

    //!< list of edge ids, which a node is incident to as Tail. For undirected
    // hypergraphs this is the same as the list of head incidence.
    std::vector<std::vector<edgeid>> nodeTailIncidence;

    //!< list of node weights
    std::vector<std::vector<nodeweight>> nodeWeights;

    // Edge related data

    //!< edgeExists[v] is false if node v has been removed from the graph
    std::vector<bool> edgeExists;

    //!< list of edge ids, which a node is incident to
    std::vector<std::vector<edgeid>> edgeIncidence;

    //!< list of edge weights
    std::vector<std::vector<edgeweight>> edgeWeights;
}

} // namespace NetworKit

#endif //
