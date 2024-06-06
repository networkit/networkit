/*
 * Graph.hpp
 *
 *  Created on: 24.05.2024
 *      Author: Fabian Brandt-Tumescheit
 */

#ifndef NETWORKIT_GRAPH_HYPERGRAPH_HPP_
#define NETWORKIT_GRAPH_HYPERGRAPH_HPP_

#include <set>
#include <vector>

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
    // bool directed;

    // Node related data

    //!< nodeExists[v] is false if node v has been removed from the graph
    std::vector<bool> nodeExists;

    //!< list of node weights
    std::vector<nodeweight> nodeWeights;

    //!< list of edge ids, which a node is incident to.
    std::vector<std::set<edgeid>> nodeIncidence;

    // Edge related data

    //!< edgeExists[v] is false if node v has been removed from the graph
    std::vector<bool> edgeExists;

    //!< list of edge weights
    std::vector<edgeweight> edgeWeights;

    //!< list of node ids, which are part of a certain hyperedge.
    std::vector<std::set<edgeid>> edgeIncidence;

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

    /**
     * Add a new node to the hypergraph and return it.
     * @return The new node.
     */
    node addNode();

    /**
     * Connect a node @a u to a list of hyperedges @a edges. If no node is given, a new
     * node is added to the hypergraph and directly connected  to a list of hyperedges
     * @a edges.
     * @return The new node.
     */
    node addNodeTo(std::vector<edgeid> edges, node u = none);

    /**
     * Remove a node @a v from the hypergraph. This removes it from every hyperedge.
     * @param v Node.
     */
    void removeNode(node v);

    /**
     * Remove a node @a v from a certain edge in  the hypergraph.
     * @param v Node.
     * @param eid Edge.
     */
    void removeNodeFrom(node v, edgeid eid);

    /* EDGE MODIFIERS */

    /**
     * Add a new edge to the hypergraph and return it.
     * @return The new edge or none, if the creation was unsuccessful
     */
    edgeid addEdge();

    /**
     * TODO: Add documentation string.
     * @param nodes Nodes to add to this edge. If a node is non-existing, this node is added if @a
     * addMissing is set to true.
     * @param ew Optional edge weight.
     * @param addMissing Adds unknown nodes if set to true. Note, that this increases the running
     * time due to additional checks.
     * @return The new edge or none, if the creation was unsuccessful
     */
    edgeid addEdge(std::vector<node> nodes, edgeweight ew = defaultEdgeWeight,
                   bool addMissing = false);

    bool removeEdge(...);

    /**
     * Return node weight of node @a u. Returns 0 if node does not
     * exist.
     *
     * @param u Node id
     * @return Node weight of @a u or 0 if edge does not exist.
     */
    nodeweight weight(node u) const;

    /**
     * Set the weight of a node. If the node does not exist,
     * it will be inserted.
     *
     * @param[in]	u	node id
     * @param[in]	weight	node weight
     */
    void setWeight(node u, nodeweight nw);
};

} // namespace NetworKit

#endif //
