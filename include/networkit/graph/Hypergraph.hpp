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

    /* BASIC ATTRIBUTES */

    // current number of nodes
    count n;

    // current number of edges
    count m;

    //!< current upper bound of node ids, z will be the id of the next node
    node z;
    //!< current upper bound of edge ids, will be the id of the next edge
    edgeid omega;

    //!< true if the hypergraph is weighted, false otherwise
    bool weighted;
    //!< true if the hypergraph is directed, false otherwise
    bool directed;

    // Node related data

    //!< nodeExists[v] is false if node v has been removed from the graph
    std::vector<bool> nodeExists;

    //!< list of node weights
    std::vector<std::vector<nodeweight>> nodeWeights;

    //!< list of edge ids, which a node is incident to as Tail. For undirected
    // hypergraphs this is the same as the list of head incidence.
    std::vector<std::vector<edgeid>> nodeTailIncidence;
    std::vector<std::vector<edgeid>> nodeHeadIncidence;

    // Edge related data

    //!< edgeExists[v] is false if node v has been removed from the graph
    std::vector<bool> edgeExists;

    //!< list of edge weights
    std::vector<std::vector<edgeweight>> edgeWeights;

    //!< list of node ids, which are part of an edge as tail. For undirected
    // hypergraphs this is the same as the list of head incidence.
    std::vector<std::vector<edgeid>> edgeTailIncidence;
    std::vector<std::vector<edgeid>> edgeHeadIncidence;

    AttributeMap<PerNode, Hypergraph> nodeAttributeMap;
    AttributeMap<PerEdge, Hypergraph> edgeAttributeMap;

public:
    /**
     * @brief Construct a new Hypergraph object
     *
     * @param n Number of nodes
     * @param m Number of edges
     * @param weighted Set the Hypergraph to be weighted
     * @param directed Set the Hypergraph to be directed
     */
    Hypergraph(count n, count m, bool weighted, bool directed);

    /* GLOBAL PROPERTIES */

    /**
     * Returns <code>true</code> if this hypergraph supports edge weights other
     * than 1.0.
     * @return <code>true</code> if this hypergraph supports edge weights other
     * than 1.0.
     */
    bool isWeighted() const noexcept { return weighted; }

    /**
     * Return @c true if this hypergraph supports directed edges.
     * @return @c true if this hypergraph supports directed edges.
     */
    bool isDirected() const noexcept { return directed; }

    /**
     * Return <code>true</code> if hypergraph contains no nodes.
     * @return <code>true</code> if hypergraph contains no nodes.
     */
    bool isEmpty() const noexcept { return !n; }

    /**
     * Return the number of nodes in the hypergraph.
     * @return The number of nodes.
     */
    count numberOfNodes() const noexcept { return n; }

    /**
     * Return the number of edges in the hypergraph.
     * @return The number of edges.
     */
    count numberOfEdges() const noexcept { return m; }

    /**
     * Get an upper bound for the node ids in the hypergraph.
     * @return An upper bound for the node ids.
     */
    index upperNodeIdBound() const noexcept { return z; }

    /**
     * Get an upper bound for the edge ids in the hypergraph.
     * @return An upper bound for the node ids.
     */
    index upperEdgeIdBound() const noexcept { return omega; }
};

} // namespace NetworKit

#endif //
