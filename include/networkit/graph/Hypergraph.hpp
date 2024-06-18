/*
 * Hypergraph.hpp
 *
 *  Created on: 24.05.2024
 *      Author: Fabian Brandt-Tumescheit
 */

#ifndef NETWORKIT_GRAPH_HYPERGRAPH_HPP_
#define NETWORKIT_GRAPH_HYPERGRAPH_HPP_

#include <cstddef>
#include <set>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/FunctionTraits.hpp>
#include <networkit/graph/Attributes.hpp>
#include <networkit/graph/EdgeIterators.hpp>
#include <networkit/graph/NeighborIterators.hpp>
#include <networkit/graph/NodeIterators.hpp>

namespace NetworKit {

/**
 * A unweighted hyperedge
 */
struct Hyperedge {
    std::set<node> nodes;

    Hyperedge() = default;

    Hyperedge(const std::vector<node> &otherNodes) {
        nodes = std::set<node>(otherNodes.begin(), otherNodes.end());
    }
};

/**
 * A weighted hyperedge
 */
struct WeightedHyperedge : Hyperedge {
    edgeweight weight;

    // Needed by cython
    WeightedHyperedge() : Hyperedge(), weight(std::numeric_limits<edgeweight>::max()) {}

    WeightedHyperedge(const std::vector<node> &otherNodes, edgeweight w)
        : Hyperedge(otherNodes), weight(w) {}
};

/**
 * @ingroup graph
 * A hypergraph (with optional weights) and parallel iterator methods.
 */
class Hypergraph final {

    /* BASIC ATTRIBUTES */

    // current number of nodes
    count numNodes;

    // current number of edges
    count numEdges;

    //!< current upper bound of node ids, z will be the id of the next node
    node maxNodeId;
    //!< current upper bound of edge ids, will be the id of the next edge
    edgeid maxEdgeId;

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
    std::vector<std::set<node>> edgeIncidence;

    AttributeMap<PerNode, Hypergraph> nodeAttributeMap;
    AttributeMap<PerEdge, Hypergraph> edgeAttributeMap;

    /**
     * Calls the given function f if it has only one argument, discards the
     * edgeweight
     */
    template <class F, void * = (void *)0>
    auto nodeLambda(F &f, node nId, nodeweight) const -> decltype(f(nId)) {
        return f(nId);
    }

    /**
     * Calls the given function f if it has only two arguments and the second
     * argument is of type edgeweight, discards the first node and the edge id
     * Note that the decltype check is not enough as edgeweight can be casted to
     * node.
     */
    template <class F,
              typename std::enable_if<
                  (Aux::FunctionTraits<F>::arity >= 1)
                  && std::is_same<nodeweight, typename Aux::FunctionTraits<F>::template arg<
                                                  1>::type>::value>::type * = (void *)0>
    auto nodeLambda(F &f, node nId, nodeweight nWeight) const -> decltype(f(nId, nWeight)) {
        return f(nId, nWeight);
    }

    template <bool hasWeights>
    inline nodeweight nodeWeightIteratorHelper(edgeid eId) const;

    /**
     * Calls the given function f if it has only one argument, discards the
     * edgeweight
     */
    template <class F, void * = (void *)0>
    auto edgeLambda(F &f, edgeid eId, edgeweight) const -> decltype(f(eId)) {
        return f(eId);
    }

    /**
     * Calls the given function f if it has only two arguments and the second
     * argument is of type edgeweight, discards the first node and the edge id
     * Note that the decltype check is not enough as edgeweight can be casted to
     * node.
     */
    template <class F,
              typename std::enable_if<
                  (Aux::FunctionTraits<F>::arity >= 1)
                  && std::is_same<edgeweight, typename Aux::FunctionTraits<F>::template arg<
                                                  1>::type>::value>::type * = (void *)0>
    auto edgeLambda(F &f, edgeid eId, edgeweight eWeight) const -> decltype(f(eId, eWeight)) {
        return f(eId, eWeight);
    }

    template <bool hasWeights>
    inline edgeweight edgeWeightIteratorHelper(edgeid eId) const;

public:
    /**
     * @brief Construct a new Hypergraph object
     *
     * @param n Number of nodes
     * @param m Number of edges
     * @param weighted Set the Hypergraph to be weighted
     */
    Hypergraph(count n = 0, count m = 0, bool weighted = false);

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
    bool isEmpty() const noexcept { return !numNodes; }

    /**
     * Return the number of nodes in the hypergraph.
     * @return The number of nodes.
     */
    count numberOfNodes() const noexcept { return numNodes; }

    /**
     * Return the number of edges in the hypergraph.
     * @return The number of edges.
     */
    count numberOfEdges() const noexcept { return numEdges; }

    /**
     * Get an upper bound for the node ids in the hypergraph.
     * @return An upper bound for the node ids.
     */
    index upperNodeIdBound() const noexcept { return maxNodeId; }

    /**
     * Get an upper bound for the edge ids in the hypergraph.
     * @return An upper bound for the node ids.
     */
    index upperEdgeIdBound() const noexcept { return maxEdgeId; }

    /* NODE PROPERTIES & MODIFIERS */

    /**
     * Add a new node to the hypergraph and return it.
     * @return The new node.
     */
    node addNode();

    /**
     * Add numberOfNewNodes new nodes.
     * @param  numberOfNewNodes Number of new nodes.
     * @return The index of the last node added.
     */
    node addNodes(count numberOfNewNodes);

    /**
     * Connect a node @a u to a list of hyperedges @a edges. If no node is given, a new
     * node is added to the hypergraph and directly connected to a list of hyperedges
     * @a edges.
     * @return The new node.
     */
    node addNodeTo(std::vector<edgeid> edges, node u = none);

    /**
     * Connect a list of nodes @a u to a hyperedge @a edgeid. If no edge is given, a new
     * hyperedge is added to the hypergraph.
     * @return The new hyperedge.
     */
    node addNodesTo(std::vector<node> nodes, edgeid eid = none);

    // TODO: add doc string
    bool hasNode(node u) const { return u < maxNodeId && nodeExists[u]; };

    // TODO: add doc string
    bool hasNode(node u, edgeid eid) const {
        return hasNode(u) && edgeIncidence[eid].find(eid) != edgeIncidence[eid].end();
    }

    /**
     * Remove a node @a u from the hypergraph. This removes it from every hyperedge.
     * @param u Node.
     */
    void removeNode(node u);

    // TODO: add docstring
    void restoreNode(node v);

    /**
     * Remove a node @a u from a certain edge in  the hypergraph.
     * @param u Node.
     * @param eid Edge.
     */
    void removeNodeFrom(node u, edgeid eid);

    /**
     * Return node weight of node @a u. Returns 0 if node does not
     * exist.
     *
     * @param u Node id
     * @return Node weight of @a u or 0 if node does not exist.
     */
    nodeweight getNodeWeight(node u) const;

    /**
     * Set the weight of a node. If the node does not exist,
     * it will be inserted.
     *
     * @param[in]	u	node id
     * @param[in]	weight	node weight
     */
    void setNodeWeight(node u, nodeweight nw);

    // TODO: add docstring
    count degree(node u) const;

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
    edgeid addEdge(const std::vector<node> &nodes, edgeweight ew = defaultEdgeWeight,
                   bool addMissing = false);

    /**
     * TODO: Add documentation string.
     */
    void removeEdge(edgeid edge);

    /**
     * Return edge weight of edge @a u. Returns 0 if edge does not
     * exist.
     *
     * @param u Edge id
     * @return Edge weight of @a u or 0 if edge does not exist.
     */
    edgeweight getEdgeWeight(edgeid u) const;

    /**
     * Set the weight of a edge. If the edge does not exist,
     * it will be inserted.
     *
     * @param[in]	u	edge id
     * @param[in]	weight	edge weight
     */
    void setEdgeWeight(edgeid u, edgeweight nw);

    // TODO: add doc string
    bool hasEdge(edgeid eid) const { return eid < maxEdgeId && edgeExists[eid]; };

    /* ITERATORS */

    /**
     * Iterate over all nodes of the graph and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L>
    void forNodes(L handle) const;

    /**
     * Iterate randomly over all nodes of the graph and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L>
    void parallelForNodes(L handle) const;

    /**
     * @brief Implementation of the for loop for all nodes, @see forNodes
     *
     * @param handle The handle that shall be executed for all nodes
     * @return void
     */
    template <bool hasWeights, typename L>
    inline void forNodesImpl(L handle) const;

    /**
     * @brief Implementation of the parallel for loop for all nodes, @see parallelForNodes
     *
     * @param handle The handle that shall be executed for all nodes
     * @return void
     */
    template <bool hasWeights, typename L>
    inline void parallelForNodesImpl(L handle) const;

    // For support of API: NetworKit::Hypergraph::NodeIterator
    using NodeIterator = NodeIteratorBase<Hypergraph>;
    // For support of API: NetworKit::Hypergraph::NodeRange
    using NodeRange = NodeRangeBase<Hypergraph>;

    /**
     * Get an iterable range over the nodes of the graph.
     *
     * @return Iterator range over the nodes of the graph.
     */
    NodeRange nodeRange() const noexcept { return NodeRange(*this); }

    /**
     * Iterate over all edges of the const graph and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameters <code>(node, node)</code>, <code>(node,
     * node, edgweight)</code>, <code>(node, node, edgeid)</code> or
     * <code>(node, node, edgeweight, edgeid)</code>.
     */
    template <typename L>
    void forEdges(L handle) const;

    /**
     * Iterate in parallel over all edges of the const graph and call @a
     * handle (lambda closure).
     *
     * @param handle Takes parameters <code>(node, node)</code> or
     * <code>(node, node, edgweight)</code>, <code>(node, node, edgeid)</code>
     * or <code>(node, node, edgeweight, edgeid)</code>.
     */
    template <typename L>
    void parallelForEdges(L handle) const;

    // For support of API: NetworKit::Hypergraph:EdgeIterator
    using EdgeIterator = EdgeTypeIterator<Hypergraph, Hyperedge>;
    // For support of API: NetworKit::Hypergraph:EdgeWeightIterator
    using EdgeWeightIterator = EdgeTypeIterator<Hypergraph, WeightedHyperedge>;
    // For support of API: NetworKit::Hypergraph:EdgeRange
    using EdgeRange = EdgeTypeRange<Hypergraph, Hyperedge>;
    // For support of API: NetworKit::Hypergraph:EdgeWeightRange
    using EdgeWeightRange = EdgeTypeRange<Hypergraph, WeightedHyperedge>;

    /**
     * @brief Implementation of the for loop for all edges, @see forEdges
     *
     * @param handle The handle that shall be executed for all edges
     * @return void
     */
    template <bool hasWeights, typename L>
    inline void forEdgeImpl(L handle) const;

    /**
     * @brief Implementation of the parallel for loop for all edges, @see parallelForEdges
     *
     * @param handle The handle that shall be executed for all edges
     * @return void
     */
    template <bool hasWeights, typename L>
    inline void parallelForEdgesImpl(L handle) const;

    /**
     * Get an iterable range over the edges of the graph.
     *
     * @return Iterator range over the edges of the graph.
     */
    EdgeRange edgeRange() const noexcept { return EdgeRange(*this); }

    /**
     * Get an iterable range over the edges of the graph and their weights.
     *
     * @return Iterator range over the edges of the graph and their weights.
     */
    EdgeWeightRange edgeWeightRange() const noexcept { return EdgeWeightRange(*this); }
};

template <typename L>
void Hypergraph::forNodes(L handle) const {
    switch (static_cast<count>(weighted)) {
    case 0: // unweighted
        forNodesImpl<false, L>(handle);
        break;

    case 1: // weighted
        forNodesImpl<true, L>(handle);
        break;
    }
}

template <typename L>
void Hypergraph::parallelForNodes(L handle) const {
    switch (static_cast<count>(weighted)) {
    case 0: // unweighted
        parallelForNodesImpl<false, L>(handle);
        break;

    case 1: // weighted
        parallelForNodesImpl<true, L>(handle);
        break;
    }
}

template <bool hasWeights, typename L>
void Hypergraph::forNodesImpl(L handle) const {
    for (node nId = 0; nId < maxNodeId; ++nId) {
        if (nodeExists[nId]) {
            nodeLambda<L>(handle, nId, nodeWeightIteratorHelper<hasWeights>(nId));
        }
    }
}

template <bool hasWeights, typename L>
void Hypergraph::parallelForNodesImpl(L handle) const {
#pragma omp parallel for
    for (omp_index nId = 0; nId < static_cast<omp_index>(maxNodeId); ++nId) {
        if (nodeExists[nId]) {
            nodeLambda<L>(handle, nId, nodeWeightIteratorHelper<hasWeights>(nId));
        }
    }
}

template <bool hasWeights, typename L>
inline void Hypergraph::forEdgeImpl(L handle) const {
    for (edgeid eId = 0; eId < maxEdgeId; ++eId) {
        if (edgeExists[eId]) {
            edgeLambda<L>(handle, eId, edgeWeightIteratorHelper<hasWeights>(eId));
        }
    }
}

template <bool hasWeights, typename L>
inline void Hypergraph::parallelForEdgesImpl(L handle) const {
#pragma omp parallel for schedule(guided)
    for (omp_index eId = 0; eId < static_cast<omp_index>(maxEdgeId); ++eId) {
        if (edgeExists[eId]) {
            edgeLambda<L>(handle, eId, edgeWeightIteratorHelper<hasWeights>(eId));
        }
    }
}

template <typename L>
void Hypergraph::forEdges(L handle) const {
    switch (static_cast<count>(weighted)) {
    case 0: // unweighted
        forEdgeImpl<false, L>(handle);
        break;

    case 1: // weighted
        forEdgeImpl<true, L>(handle);
        break;
    }
}

template <typename L>
void Hypergraph::parallelForEdges(L handle) const {
    switch (static_cast<count>(weighted)) {
    case 0: // unweighted
        parallelForEdgesImpl<false, L>(handle);
        break;

    case 1: // weighted
        parallelForEdgesImpl<true, L>(handle);
        break;
    }
}

/* HELPERS */

// implementation for weighted == true
template <bool hasWeights>
inline edgeweight Hypergraph::edgeWeightIteratorHelper(edgeid eId) const {
    return edgeWeights[eId];
}

// implementation for weighted == false
template <>
inline edgeweight Hypergraph::edgeWeightIteratorHelper<false>(edgeid) const {
    return defaultEdgeWeight;
}

// implementation for weighted == true
template <bool hasWeights>
inline nodeweight Hypergraph::nodeWeightIteratorHelper(node nId) const {
    return nodeWeights[nId];
}

// implementation for weighted == false
template <>
inline nodeweight Hypergraph::nodeWeightIteratorHelper<false>(node) const {
    return defaultEdgeWeight;
}

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_HYPERGRAPH_HPP_
