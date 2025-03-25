#ifndef NETWORKIT_GRAPH_GRAPH_HPP_
#define NETWORKIT_GRAPH_GRAPH_HPP_

#pragma message( \
    "warning: Graph.hpp is deprecated. Include a graph implementation like AdjListGraph.hpp instead.")

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/ArrayTools.hpp>
#include <networkit/auxiliary/FunctionTraits.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/Attributes.hpp>
#include <networkit/graph/Edge.hpp>
#include <networkit/graph/EdgeIterators.hpp>
#include <networkit/graph/EdgeUtils.hpp>
#include <networkit/graph/NeighborIterators.hpp>
#include <networkit/graph/NodeIterators.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

// forward declaration to randomization/CurveballImpl.hpp
namespace CurveballDetails {
class CurveballMaterialization;
}

/**
 * @ingroup graph
 * A graph (with optional weights) and parallel iterator methods.
 */
template <class NodeT, class EdgeWeightT>
class AdjListGraph final {
    static_assert(std::is_integral_v<NodeT>, "NodeT must be an integer type.");
    static_assert(std::is_arithmetic_v<EdgeWeightT>, "EdgeWeightT must be a numeric type.");
    // graph attributes
    //!< current number of nodes
    count n;
    //!< current number of edges
    count m;

    //!< current number of self loops, edges which have the same origin and
    //!< target
    count storedNumberOfSelfLoops;

    //!< current upper bound of node ids, z will be the id of the next node
    NodeT z;
    //!< current upper bound of edge ids, will be the id of the next edge
    edgeid omega;
    //!< current time step
    count t;

    //!< true if the graph is weighted, false otherwise
    bool weighted;
    //!< true if the graph is directed, false otherwise
    bool directed;
    //!< true if edge ids have been assigned
    bool edgesIndexed;

    //!< true if edge removals should maintain compact edge ids
    bool maintainCompactEdges = false;
    //!< true if edge removals should maintain sorted edge ids
    bool maintainSortedEdges = false;

    //!< saves the ID of the most recently removed edge (if exists)
    edgeid deletedID;

    // per node data
    //!< exists[v] is true if node v has not been removed from the graph
    std::vector<bool> exists;

    //!< only used for directed graphs, inEdges[v] contains all nodes u that
    //!< have an edge (u, v)
    std::vector<std::vector<NodeT>> inEdges;
    //!< (outgoing) edges, for each edge (u, v) v is saved in outEdges[u] and
    //!< for undirected also u in outEdges[v]
    std::vector<std::vector<NodeT>> outEdges;

    //!< only used for directed graphs, same schema as inEdges
    std::vector<std::vector<EdgeWeightT>> inEdgeWeights;
    //!< same schema (and same order!) as outEdges
    std::vector<std::vector<EdgeWeightT>> outEdgeWeights;

    //!< only used for directed graphs, same schema as inEdges
    std::vector<std::vector<edgeid>> inEdgeIds;
    //!< same schema (and same order!) as outEdges
    std::vector<std::vector<edgeid>> outEdgeIds;

    static constexpr NodeT nullNodeId = NullNodeId<NodeT>;

private:
    AttributeMap<PerNode, AdjListGraph> nodeAttributeMap;
    AttributeMap<PerEdge, AdjListGraph> edgeAttributeMap;

public:
    auto& nodeAttributes() noexcept { return nodeAttributeMap; }
    auto const& nodeAttributes() const noexcept { return nodeAttributeMap; }
    auto& edgeAttributes() noexcept { return edgeAttributeMap; }
    auto const& edgeAttributes() const noexcept { return edgeAttributeMap; }

    // wrap up some typed attributes for the cython interface:
    //

    auto attachNodeIntAttribute(std::string const& name) {
        nodeAttributes().theGraph = this;
        return nodeAttributes().template attach<int>(name);
    }

    auto attachEdgeIntAttribute(std::string const& name) {
        edgeAttributes().theGraph = this;
        return edgeAttributes().template attach<int>(name);
    }

    auto attachNodeDoubleAttribute(std::string const& name) {
        nodeAttributes().theGraph = this;
        return nodeAttributes().template attach<double>(name);
    }

    auto attachEdgeDoubleAttribute(std::string const& name) {
        edgeAttributes().theGraph = this;
        return edgeAttributes().template attach<double>(name);
    }

    auto attachNodeStringAttribute(std::string const& name) {
        nodeAttributes().theGraph = this;
        return nodeAttributes().template attach<std::string>(name);
    }

    auto attachEdgeStringAttribute(std::string const& name) {
        edgeAttributes().theGraph = this;
        return edgeAttributes().template attach<std::string>(name);
    }

    auto getNodeIntAttribute(std::string const& name) {
        nodeAttributes().theGraph = this;
        return nodeAttributes().template get<int>(name);
    }

    auto getEdgeIntAttribute(std::string const& name) {
        edgeAttributes().theGraph = this;
        return edgeAttributes().template get<int>(name);
    }

    auto getNodeDoubleAttribute(std::string const& name) {
        nodeAttributes().theGraph = this;
        return nodeAttributes().template get<double>(name);
    }

    auto getEdgeDoubleAttribute(std::string const& name) {
        edgeAttributes().theGraph = this;
        return edgeAttributes().template get<double>(name);
    }

    auto getNodeStringAttribute(std::string const& name) {
        nodeAttributes().theGraph = this;
        return nodeAttributes().template get<std::string>(name);
    }

    auto getEdgeStringAttribute(std::string const& name) {
        edgeAttributes().theGraph = this;
        return edgeAttributes().template get<std::string>(name);
    }

    void detachNodeAttribute(std::string const& name) {
        nodeAttributes().theGraph = this;
        nodeAttributes().detach(name);
    }

    void detachEdgeAttribute(std::string const& name) {
        edgeAttributes().theGraph = this;
        edgeAttributes().detach(name);
    }

    using NodeIntAttribute = Attribute<PerNode, AdjListGraph, int, false>;
    using NodeDoubleAttribute = Attribute<PerNode, AdjListGraph, double, false>;
    using NodeStringAttribute = Attribute<PerNode, AdjListGraph, std::string, false>;

    using EdgeIntAttribute = Attribute<PerEdge, AdjListGraph, int, false>;
    using EdgeDoubleAttribute = Attribute<PerEdge, AdjListGraph, double, false>;
    using EdgeStringAttribute = Attribute<PerEdge, AdjListGraph, std::string, false>;

private:
    static constexpr EdgeWeightT defaultEdgeWeightT = EdgeWeightT{1};
    /**
     * Returns the index of node u in the array of incoming edges of node v.
     * (for directed graphs inEdges is searched, while for indirected outEdges
     * is searched, which gives the same result as indexInOutEdgeArray).
     */
    index indexInInEdgeArray(NodeT v, NodeT u) const;

    /**
     * Returns the index of node v in the array of outgoing edges of node u.
     */
    index indexInOutEdgeArray(NodeT u, NodeT v) const;

    /**
     * Computes the weighted in/out degree of node @a u.
     *
     * @param u Node.
     * @param inDegree whether to compute the in degree or the out degree.
     * @param countSelfLoopsTwice If set to true, self-loops will be counted twice.
     *
     * @return Weighted in/out degree of node @a u.
     */
    EdgeWeightT computeWeightedDegree(NodeT u, bool inDegree = false,
                                      bool countSelfLoopsTwice = false) const;

    /**
     * Returns the edge weight of the outgoing edge of index i in the outgoing
     * edges of node u
     * @param u The node
     * @param i The index
     * @return The weight of the outgoing edge or defaultEdgeWeightT if the graph
     * is unweighted
     */
    template <bool hasWeights>
    inline EdgeWeightT getOutEdgeWeight(NodeT u, index i) const {
        if constexpr (hasWeights)
            return outEdgeWeights[u][i];
        else
            return defaultEdgeWeightT;
    }

    /**
     * Returns the edge weight of the incoming edge of index i in the incoming
     * edges of node u
     *
     * @param u The node
     * @param i The index in the incoming edge array
     * @return The weight of the incoming edge
     */
    template <bool hasWeights>
    inline EdgeWeightT getInEdgeWeight(NodeT u, index i) const {
        if constexpr (hasWeights)
            return inEdgeWeights[u][i];
        else
            return defaultEdgeWeightT;
    }

    /**
     * Returns the edge id of the edge of index i in the outgoing edges of node
     * u
     *
     * @param u The node
     * @param i The index in the outgoing edges
     * @return The edge id
     */
    template <bool graphHasEdgeIds>
    inline edgeid getOutEdgeId(NodeT u, index i) const {
        if constexpr (graphHasEdgeIds)
            return outEdgeIds[u][i];
        else
            return nullEdgeId;
    }

    /**
     * Returns the edge id of the edge of index i in the incoming edges of node
     * u
     *
     * @param u The node
     * @param i The index in the incoming edges of u
     * @return The edge id
     */
    template <bool graphHasEdgeIds>
    inline edgeid getInEdgeId(NodeT u, index i) const {
        if constexpr (graphHasEdgeIds)
            return inEdgeIds[u][i];
        else
            return nullEdgeId;
    }

    /**
     * @brief Returns if the edge (u, v) shall be used in the iteration of all
     * edgesIndexed
     *
     * @param u The source node of the edge
     * @param v The target node of the edge
     * @return If the node shall be used, i.e. if v is not nullNodeId and in the
     * undirected case if u >= v
     */
    template <bool graphIsDirected>
    inline bool useEdgeInIteration(NodeT u, NodeT v) const {
        if constexpr (graphIsDirected)
            return true;
        else
            return u >= v;
    }

    /**
     * @brief Implementation of the for loop for outgoing edges of u
     *
     * Note: If all (valid) outgoing edges shall be considered, graphIsDirected
     * needs to be set to true
     *
     * @param u The node
     * @param handle The handle that shall be executed for each edge
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
    inline void forOutEdgesOfImpl(NodeT u, L handle) const;

    /**
     * @brief Implementation of the for loop for incoming edges of u
     *
     * For undirected graphs, this is the same as forOutEdgesOfImpl but u and v
     * are changed in the handle
     *
     * @param u The node
     * @param handle The handle that shall be executed for each edge
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
    inline void forInEdgesOfImpl(NodeT u, L handle) const;

    /**
     * @brief Implementation of the for loop for all edges, @see forEdges
     *
     * @param handle The handle that shall be executed for all edges
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
    inline void forEdgeImpl(L handle) const;

    /**
     * @brief Parallel implementation of the for loop for all edges, @see
     * parallelForEdges
     *
     * @param handle The handle that shall be executed for all edges
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
    inline void parallelForEdgesImpl(L handle) const;

    /**
     * @brief Summation variant of the parallel for loop for all edges, @see
     * parallelSumForEdges
     *
     * @param handle The handle that shall be executed for all edges
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
    inline double parallelSumForEdgesImpl(L handle) const;

    /* Compile-time callable dispatcher.
     *
     * Aux::FunctionTraits lets us inspect a callable's arity and parameter types
     * so we can invoke it with the best matching edge callback signature.
     * SafeArg prevents out-of-bounds access to Traits::arg<N> by returning void
     * when the callable has fewer than N + 1 parameters.
     * ArgType normalizes the inspected type with std::decay_t so comparisons work
     * for references, const qualifiers, and other call-compatible forms.
     *
     * The overload selection below is ordered from most specific to least specific.
     * This matters because several lambda signatures overlap, especially when node
     * and EdgeWeightT are the same type.
     *
     * If no branch matches, the static_assert produces a clear compile-time error
     * instead of a long substitution failure.
     */

    // Return void when the callable does not have an Nth parameter.
    template <typename TraitsT, size_t N, typename Enable = void>
    struct SafeArg {
        using type = void;
    };

    // Enable this specialization only when N is in range.
    template <typename TraitsT, size_t N>
    struct SafeArg<TraitsT, N, typename std::enable_if<(N < TraitsT::arity)>::type> {
        using type = typename TraitsT::template arg<N>::type;
    };

    template <typename F, size_t N>
    using ArgType = std::decay_t<typename SafeArg<Aux::FunctionTraits<F>, N>::type>;

    template <class F>
    decltype(auto) edgeLambda(F& f, NodeT u, NodeT v, EdgeWeightT ew, edgeid id) const {
        using Traits = Aux::FunctionTraits<F>;
        constexpr size_t arity = Traits::arity;

        // Most specific signatures first to avoid ambiguity between overlapping
        // callable shapes and implicit conversions.
        if constexpr (arity >= 3 && std::is_same_v<ArgType<F, 2>, EdgeWeightT>
                      && std::is_same_v<ArgType<F, 3>, edgeid>) {
            return f(u, v, ew, id);
        } else if constexpr (arity >= 2 && std::is_same_v<ArgType<F, 2>, edgeid>
                             && std::is_same_v<ArgType<F, 1>, NodeT>) {
            return f(u, v, id);
        } else if constexpr (arity >= 2 && std::is_same_v<ArgType<F, 2>, EdgeWeightT>) {
            return f(u, v, ew);
        } else if constexpr (arity >= 1 && std::is_same_v<ArgType<F, 1>, EdgeWeightT>
                             && !std::is_same_v<NodeT, EdgeWeightT>) {
            // Treat this as a weight-taking callback only when node ids and weights
            // are distinct types.
            return f(v, ew);
        } else if constexpr (arity >= 1 && std::is_same_v<ArgType<F, 1>, NodeT>) {
            return f(u, v);
        } else if constexpr (arity >= 1 && std::is_same_v<ArgType<F, 1>, EdgeWeightT>) {
            return f(v, ew);
        } else if constexpr (arity >= 1 && std::is_same_v<ArgType<F, 0>, NodeT>) {
            return f(v);
        } else {
            static_assert(!std::is_same_v<F, F>,
                          "Your lambda does not support the required parameters or the parameters "
                          "have the wrong type.");
        }
    }

    /**
     * Calls the given BFS handle with distance parameter
     */
    template <class F>
    auto callBFSHandle(F& f, NodeT u, count dist) const -> decltype(f(u, dist)) {
        return f(u, dist);
    }

    /**
     * Calls the given BFS handle without distance parameter
     */
    template <class F>
    auto callBFSHandle(F& f, NodeT u, count) const -> decltype(f(u)) {
        return f(u);
    }

public:
    // For support of API: NetworKit::AdjListGraph::NodeIterator
    using NodeIterator = NodeIteratorBase<AdjListGraph, NodeT, EdgeWeightT>;
    // For support of API: NetworKit::AdjListGraph::NodeRange
    using NodeRange = NodeRangeBase<AdjListGraph, NodeT, EdgeWeightT>;

    // For support of API: NetworKit::AdjListGraph:EdgeIterator
    using EdgeIterator = EdgeWeightTIterator<AdjListGraph, NodeT, EdgeWeightT, EdgeT<NodeT>>;
    // For support of API: NetworKit::AdjListGraph:EdgeWeightIterator
    using EdgeWeightIterator =
        EdgeWeightTIterator<AdjListGraph, NodeT, EdgeWeightT, WeightedEdgeT<NodeT, EdgeWeightT>>;
    // For support of API: NetworKit::AdjListGraph:EdgeRange
    using EdgeRange = EdgeWeightTRange<AdjListGraph, NodeT, EdgeWeightT, EdgeT<NodeT>>;
    // For support of API: NetworKit::AdjListGraph:EdgeWeightRange
    using EdgeWeightRange =
        EdgeWeightTRange<AdjListGraph, NodeT, EdgeWeightT, WeightedEdgeT<NodeT, EdgeWeightT>>;

    // For support of API: NetworKit::AdjListGraph::NeighborIterator;
    using NeighborIterator = NeighborIteratorBase<std::vector<NodeT>>;
    // For support of API: NetworKit::AdjListGraph::NeighborIterator;
    using NeighborWeightIterator =
        NeighborWeightIteratorBase<std::vector<NodeT>, std::vector<EdgeWeightT>>;

    /**
     * Wrapper class to iterate over a range of the neighbors of a node within
     * a for loop.
     */
    template <bool InEdges = false>
    class NeighborRange {
        AdjListGraph const* G;
        NodeT u{nullNodeId};

    public:
        NeighborRange(AdjListGraph const& G, NodeT u) : G(&G), u(u) { assert(G.hasNode(u)); };

        NeighborRange() : G(nullptr){};

        NeighborIterator begin() const {
            assert(G);
            return InEdges ? NeighborIterator(G->inEdges[u].begin())
                           : NeighborIterator(G->outEdges[u].begin());
        }

        NeighborIterator end() const {
            assert(G);
            return InEdges ? NeighborIterator(G->inEdges[u].end())
                           : NeighborIterator(G->outEdges[u].end());
        }
    };

    using OutNeighborRange = NeighborRange<false>;

    using InNeighborRange = NeighborRange<true>;
    /**
     * Wrapper class to iterate over a range of the neighbors of a node
     * including the edge weights within a for loop.
     * Values are std::pair<NodeT, edgeweight>.
     */
    template <bool InEdges = false>
    class NeighborWeightRange {

        AdjListGraph const* G;
        NodeT u{nullNodeId};

    public:
        NeighborWeightRange(AdjListGraph const& G, NodeT u) : G(&G), u(u) { assert(G.hasNode(u)); };

        NeighborWeightRange() : G(nullptr){};

        NeighborWeightIterator begin() const {
            assert(G);
            return InEdges
                       ? NeighborWeightIterator(G->inEdges[u].begin(), G->inEdgeWeights[u].begin())
                       : NeighborWeightIterator(G->outEdges[u].begin(),
                                                G->outEdgeWeights[u].begin());
        }

        NeighborWeightIterator end() const {
            assert(G);
            return InEdges
                       ? NeighborWeightIterator(G->inEdges[u].end(), G->inEdgeWeights[u].end())
                       : NeighborWeightIterator(G->outEdges[u].end(), G->outEdgeWeights[u].end());
        }
    };

    using OutNeighborWeightRange = NeighborWeightRange<false>;

    using InNeighborWeightRange = NeighborWeightRange<true>;

    /**
     * Create a graph of @a n nodes. The graph has assignable edge weights if @a
     * weighted is set to <code>true</code>. If @a weighted is set to
     * <code>false</code> each edge has edge weight 1.0 and any other weight
     * assignment will be ignored.
     * @param n Number of nodes.
     * @param weighted If set to <code>true</code>, the graph has edge weights.
     * @param directed If set to @c true, the graph will be directed.
     */
    AdjListGraph(count n = 0, bool weighted = false, bool directed = false,
                 bool edgesIndexed = false);

    template <class EdgeMerger = std::plus<EdgeWeightT>>
    AdjListGraph(AdjListGraph const& G, bool weighted, bool directed, bool edgesIndexed = false,
                 EdgeMerger edgeMerger = std::plus<EdgeWeightT>())
        : n(G.n), m(G.m), storedNumberOfSelfLoops(G.storedNumberOfSelfLoops), z(G.z),
          omega(edgesIndexed ? G.omega : 0), t(G.t), weighted(weighted), directed(directed),
          edgesIndexed(edgesIndexed), // edges are not indexed by default
          exists(G.exists),

          // let the following be empty for the start, we fill them later
          inEdges(0), outEdges(0), inEdgeWeights(0), outEdgeWeights(0), inEdgeIds(0), outEdgeIds(0),

          // copy node attribute map
          nodeAttributeMap(G.nodeAttributeMap, this),
          // fill this later
          edgeAttributeMap(this) {

        if (G.isDirected() == directed) {
            // G.inEdges might be empty (if G is undirected), but
            // that's fine
            inEdges = G.inEdges;
            outEdges = G.outEdges;
            edgeAttributeMap = AttributeMap(G.edgeAttributeMap, this);

            // copy weights if needed
            if (weighted) {
                if (G.isWeighted()) {
                    // just copy from G, again either both graphs are directed or both are
                    // undirected
                    inEdgeWeights = G.inEdgeWeights;
                    outEdgeWeights = G.outEdgeWeights;
                } else {
                    // G has no weights, set defaultEdgeWeightT for all edges
                    if (directed) {
                        inEdgeWeights.resize(z);
                        for (NodeT u = 0; u < z; ++u) {
                            inEdgeWeights[u].resize(G.inEdges[u].size(), defaultEdgeWeightT);
                        }
                    }

                    outEdgeWeights.resize(z);
                    for (NodeT u = 0; u < z; ++u) {
                        outEdgeWeights[u].resize(outEdges[u].size(), defaultEdgeWeightT);
                    }
                }
            }
            if (G.hasEdgeIds() && edgesIndexed) {
                inEdgeIds = G.inEdgeIds;
                outEdgeIds = G.outEdgeIds;
            }
        } else if (G.isDirected()) {
            // G is directed, but we want an undirected graph
            // so we need to combine the out and in stuff for every node

            // for edge attributes, it is not well defined how they should be comined - we will skip
            // that step and warn the user
            WARN("Edge attributes are not preserved when converting from directed to undirected "
                 "graphs. The resulting graph will have empty edge attributes.");

            outEdges.resize(z);
            if (weighted)
                outEdgeWeights.resize(z);
            if (G.hasEdgeIds() && edgesIndexed)
                outEdgeIds.resize(z);
            G.balancedParallelForNodes([&](NodeT u) {
                // copy both out and in edges into our new outEdges
                outEdges[u].reserve(G.outEdges[u].size() + G.inEdges[u].size());
                outEdges[u].insert(outEdges[u].end(), G.outEdges[u].begin(), G.outEdges[u].end());
                if (weighted) {
                    if (G.isWeighted()) {
                        // same for weights
                        outEdgeWeights[u].reserve(G.outEdgeWeights[u].size()
                                                  + G.inEdgeWeights[u].size());
                        outEdgeWeights[u].insert(outEdgeWeights[u].end(),
                                                 G.outEdgeWeights[u].begin(),
                                                 G.outEdgeWeights[u].end());
                    } else {
                        // we are undirected, so no need to write anything into inEdgeWeights
                        outEdgeWeights[u].resize(outEdges[u].size(), defaultEdgeWeightT);
                    }
                }
                if (G.hasEdgeIds() && edgesIndexed) {
                    // copy both out and in edges ids into our new outEdgesIds
                    outEdgeIds[u].reserve(G.outEdgeIds[u].size() + G.inEdgeIds[u].size());
                    outEdgeIds[u].insert(outEdgeIds[u].end(), G.outEdgeIds[u].begin(),
                                         G.outEdgeIds[u].end());
                }
            });
            G.balancedParallelForNodes([&](NodeT u) {
                // this is necessary to avoid multi edges, because both u -> v and v -> u can exist
                // in G
                count edgeSurplus = 0;
                for (count i = 0; i < G.inEdges[u].size(); ++i) {
                    NodeT v = G.inEdges[u][i];
                    bool alreadyPresent = false;
                    for (count j = 0; j < G.outEdges[u].size(); ++j) {
                        if (v != G.outEdges[u][j])
                            continue; // the edge already exists as an out edge
                        alreadyPresent = true;
                        if (u != v) {
                            ++edgeSurplus;
                            if (weighted) // we need combine those edges weights when making it a
                                          // single edge
                                outEdgeWeights[u][j] =
                                    G.isWeighted()
                                        ? edgeMerger(G.inEdgeWeights[u][i], G.outEdgeWeights[u][j])
                                        : edgeMerger(defaultEdgeWeightT, defaultEdgeWeightT);
                            if (G.hasEdgeIds() && edgesIndexed)
                                outEdgeIds[u][j] = std::min(G.inEdgeIds[u][i], G.outEdgeIds[u][j]);
                        }
                        break;
                    }
                    if (!alreadyPresent) { // an equivalent out edge wasn't present so we add it
                        outEdges[u].push_back(v);
                        if (weighted)
                            outEdgeWeights[u].push_back(G.isWeighted() ? G.inEdgeWeights[u][i]
                                                                       : defaultEdgeWeightT);
                        if (G.hasEdgeIds() && edgesIndexed)
                            outEdgeIds[u].push_back(G.inEdgeIds[u][i]);
                    }
                }
#pragma omp atomic
                m -= edgeSurplus;
            });
        } else {
            // G is not directed, but this copy should be
            // generally we can can copy G.out stuff into our in stuff

            // for edge attributes, we currently do not have a way to iterate them for duplication
            WARN("Edge attributes are currently not preserved when converting from undirected to "
                 "directed graphs.");

            inEdges = G.outEdges;
            outEdges = G.outEdges;
            if (weighted) {
                if (G.isWeighted()) {
                    inEdgeWeights = G.outEdgeWeights;
                    outEdgeWeights = G.outEdgeWeights;
                } else {
                    // initialize both inEdgeWeights and outEdgeWeights with the
                    // defaultEdgeWeightT
                    inEdgeWeights.resize(z);
                    for (NodeT u = 0; u < z; ++u) {
                        inEdgeWeights[u].resize(inEdges[u].size(), defaultEdgeWeightT);
                    }
                    outEdgeWeights.resize(z);
                    for (NodeT u = 0; u < z; ++u) {
                        outEdgeWeights[u].resize(outEdges[u].size(), defaultEdgeWeightT);
                    }
                }
            }
            if (G.hasEdgeIds() && edgesIndexed) {
                inEdgeIds = G.outEdgeIds;
                outEdgeIds = G.outEdgeIds;
            }
        }

        if (!G.edgesIndexed && edgesIndexed)
            indexEdges();
    }

    /**
     * Generate a weighted graph from a list of edges. (Useful for small
     * graphs in unit tests that you do not want to read from a file.)
     *
     * @param[in] edges list of weighted edges
     */
    AdjListGraph(std::initializer_list<WeightedEdgeT<NodeT, EdgeWeightT>> edges);

    /**
     * Create a graph as copy of @a other.
     * @param other The graph to copy.
     */
    AdjListGraph(AdjListGraph const& other)
        : n(other.n), m(other.m), storedNumberOfSelfLoops(other.storedNumberOfSelfLoops),
          z(other.z), omega(other.omega), t(other.t), weighted(other.weighted),
          directed(other.directed), edgesIndexed(other.edgesIndexed), deletedID(other.deletedID),
          exists(other.exists), inEdges(other.inEdges), outEdges(other.outEdges),
          inEdgeWeights(other.inEdgeWeights), outEdgeWeights(other.outEdgeWeights),
          inEdgeIds(other.inEdgeIds), outEdgeIds(other.outEdgeIds),
          // call special constructors to copy attribute maps
          nodeAttributeMap(other.nodeAttributeMap, this),
          edgeAttributeMap(other.edgeAttributeMap, this){};

    /** move constructor */
    AdjListGraph(AdjListGraph&& other) noexcept
        : n(other.n), m(other.m), storedNumberOfSelfLoops(other.storedNumberOfSelfLoops),
          z(other.z), omega(other.omega), t(other.t), weighted(other.weighted),
          directed(other.directed), edgesIndexed(other.edgesIndexed), deletedID(other.deletedID),
          exists(std::move(other.exists)), inEdges(std::move(other.inEdges)),
          outEdges(std::move(other.outEdges)), inEdgeWeights(std::move(other.inEdgeWeights)),
          outEdgeWeights(std::move(other.outEdgeWeights)), inEdgeIds(std::move(other.inEdgeIds)),
          outEdgeIds(std::move(other.outEdgeIds)),
          nodeAttributeMap(std::move(other.nodeAttributeMap)),
          edgeAttributeMap(std::move(other.edgeAttributeMap)) {
        // attributes: set graph pointer to this new graph
        nodeAttributeMap.theGraph = this;
        edgeAttributeMap.theGraph = this;
    };

    /** Default destructor */
    ~AdjListGraph() = default;

    /** move assignment operator */
    AdjListGraph& operator=(AdjListGraph&& other) noexcept {
        std::swap(n, other.n);
        std::swap(m, other.m);
        std::swap(storedNumberOfSelfLoops, other.storedNumberOfSelfLoops);
        std::swap(z, other.z);
        std::swap(omega, other.omega);
        std::swap(t, other.t);
        std::swap(weighted, other.weighted);
        std::swap(directed, other.directed);
        std::swap(edgesIndexed, other.edgesIndexed);
        std::swap(exists, other.exists);
        std::swap(inEdges, other.inEdges);
        std::swap(outEdges, other.outEdges);
        std::swap(inEdgeWeights, other.inEdgeWeights);
        std::swap(outEdgeWeights, other.outEdgeWeights);
        std::swap(inEdgeIds, other.inEdgeIds);
        std::swap(outEdgeIds, other.outEdgeIds);
        std::swap(deletedID, other.deletedID);

        // attributes: set graph pointer to this new graph
        std::swap(nodeAttributeMap, other.nodeAttributeMap);
        std::swap(edgeAttributeMap, other.edgeAttributeMap);
        nodeAttributeMap.theGraph = this;
        edgeAttributeMap.theGraph = this;

        return *this;
    };

    /** copy assignment operator */
    AdjListGraph& operator=(AdjListGraph const& other) {
        n = other.n;
        m = other.m;
        storedNumberOfSelfLoops = other.storedNumberOfSelfLoops;
        z = other.z;
        omega = other.omega;
        t = other.t;
        weighted = other.weighted;
        directed = other.directed;
        edgesIndexed = other.edgesIndexed;
        exists = other.exists;
        inEdges = other.inEdges;
        outEdges = other.outEdges;
        inEdgeWeights = other.inEdgeWeights;
        outEdgeWeights = other.outEdgeWeights;
        inEdgeIds = other.inEdgeIds;
        outEdgeIds = other.outEdgeIds;
        deletedID = other.deletedID;

        // call special constructors to copy attribute maps
        nodeAttributeMap = AttributeMap(other.nodeAttributeMap, this);
        edgeAttributeMap = AttributeMap(other.edgeAttributeMap, this);

        return *this;
    };

    /**
     * Reserves memory in the node's edge containers for undirected graphs.
     *
     * @param u the node memory should be reserved for
     * @param size the amount of memory to reserve
     *
     * This function is thread-safe if called from different
     * threads on different nodes.
     */
    void preallocateUndirected(NodeT u, size_t size);

    /**
     * Reserves memory in the node's edge containers for directed graphs.
     *
     * @param u the node memory should be reserved for
     * @param inSize the amount of memory to reserve for in edges
     * @param outSize the amount of memory to reserve for out edges
     *
     * This function is thread-safe if called from different
     * threads on different nodes.
     */
    void preallocateDirected(NodeT u, size_t outSize, size_t inSize);

    /**
     * Reserves memory in the node's edge containers for directed graphs.
     *
     * @param u the node memory should be reserved for
     * @param outSize the amount of memory to reserve for out edges
     *
     * This function is thread-safe if called from different
     * threads on different nodes.
     */
    void preallocateDirectedOutEdges(NodeT u, size_t outSize);

    /**
     * Reserves memory in the node's edge containers for directed graphs.
     *
     * @param u the node memory should be reserved for
     * @param inSize the amount of memory to reserve for in edges
     *
     * This function is thread-safe if called from different
     * threads on different nodes.
     */
    void preallocateDirectedInEdges(NodeT u, size_t inSize);

    /** EDGE IDS **/

    /**
     * Initially assign integer edge identifiers.
     *
     * @param force Force re-indexing of edges even if they have already been
     * indexed
     */
    void indexEdges(bool force = false);

    /**
     * Checks if edges have been indexed
     *
     * @return bool if edges have been indexed
     */
    bool hasEdgeIds() const noexcept { return edgesIndexed; }

    /**
     * Get the id of the given edge.
     */
    edgeid edgeId(NodeT u, NodeT v) const;

    /**
     * Get the Edge (u,v) of the given id. (inverse to edgeId)
     * @note Time complexity of this function is O(n).
     */
    std::pair<NodeT, NodeT> edgeById(index id) const {
        std::pair<NodeT, NodeT> result{nullNodeId, nullNodeId};
        bool found = false;

        forNodesWhile([&] { return !found; },
                      [&](NodeT u) {
                          forNeighborsOf(u, [&](NodeT v) {
                              if (!this->isDirected() && v < u)
                                  return;
                              auto uvId = edgeId(u, v);
                              if (uvId == id) {
                                  found = true;
                                  result = std::make_pair(u, v);
                              }
                          });
                      });

        return result;
    }

    /**
     * Get an upper bound for the edge ids in the graph.
     * @return An upper bound for the edge ids.
     */
    index upperEdgeIdBound() const noexcept { return omega; }

    /** GRAPH INFORMATION **/

    /**
     * Try to save some memory by shrinking internal data structures of the
     * graph. Only run this once you finished editing the graph. Otherwise it
     * will cause unnecessary reallocation of memory.
     */
    void shrinkToFit();

    /**
     * DEPRECATED: this function will no longer be supported in later releases.
     * Compacts the adjacency arrays by re-using no longer needed slots from
     * deleted edges.
     */
    void TLX_DEPRECATED(compactEdges());

    /**
     * Sorts the outgoing neighbors of a given node according to a user-defined comparison function.
     *
     * @param u The node whose outgoing neighbors will be sorted.
     * @param lambda A binary predicate used to compare two neighbors. The predicate should
     *               take two nodes as arguments and return true if the first node should
     *               precede the second in the sorted order.
     */
    template <typename Lambda>
    void sortNeighbors(NodeT u, Lambda lambda);

    /**
     * Sorts the adjacency arrays by node id. While the running time is linear
     * this temporarily duplicates the memory.
     */
    void sortEdges();

    /**
     * Sorts the adjacency arrays by a custom criterion.
     *
     * @param lambda Lambda function used to sort the edges. It takes two WeightedEdge
     * e1 and e2 as input parameters, returns true if e1 < e2, false otherwise.
     */
    template <class Lambda>
    void sortEdges(Lambda lambda);

    /**
     * Set edge count of the graph to edges.
     * @param edges the edge count of a graph
     */
    void setEdgeCount(Unsafe, count edges) { m = edges; }

    /**
     * Set upper bound of edge count.
     *
     * @param newBound New upper edge id bound.
     */
    void setUpperEdgeIdBound(Unsafe, edgeid newBound) { omega = newBound; }

    /**
     * Set the number of self-loops.
     *
     * @param loops New number of self-loops.
     */
    void setNumberOfSelfLoops(Unsafe, count loops) { storedNumberOfSelfLoops = loops; }

    /* NODE MODIFIERS */

    /**
     * Add a new node to the graph and return it.
     * @return The new node.
     */
    NodeT addNode();

    /**
     * Add numberOfNewNodes new nodes.
     * @param  numberOfNewNodes Number of new nodes.
     * @return The index of the last node added.
     */
    NodeT addNodes(count numberOfNewNodes);

    /**
     * Remove a node @a v and all incident edges from the graph.
     *
     * Incoming as well as outgoing edges will be removed.
     *
     * @param v Node.
     */
    void removeNode(NodeT v);

    /**
     * Removes out-going edges from node @u. If the graph is weighted and/or has edge ids, weights
     * and/or edge ids will also be removed.
     *
     * @param u Node.
     */
    void removePartialOutEdges(Unsafe, NodeT u) {
        assert(hasNode(u));
        outEdges[u].clear();
        if (isWeighted()) {
            outEdgeWeights[u].clear();
        }
        if (hasEdgeIds()) {
            outEdgeIds[u].clear();
        }
    }

    /**
     * Removes in-going edges to node @u. If the graph is weighted and/or has edge ids, weights
     * and/or edge ids will also be removed.
     *
     * @param u Node.
     */
    void removePartialInEdges(Unsafe, NodeT u) {
        assert(hasNode(u));
        inEdges[u].clear();
        if (isWeighted()) {
            inEdgeWeights[u].clear();
        }
        if (hasEdgeIds()) {
            inEdgeIds[u].clear();
        }
    }

    /**
     * Check if node @a v exists in the graph.
     *
     * @param v Node.
     * @return @c true if @a v exists, @c false otherwise.
     */

    bool hasNode(NodeT v) const noexcept { return (v < z) && this->exists[v]; }

    /**
     * Restores a previously deleted node @a v with its previous id in the
     * graph.
     *
     * @param v Node.
     *
     */

    void restoreNode(NodeT v);

    /** NODE PROPERTIES **/
    /**
     * Returns the number of outgoing neighbors of @a v.
     *
     * @param v Node.
     * @return The number of outgoing neighbors.
     * @note The existence of the node is not checked. Calling this function with a non-existing
     * node results in a segmentation fault. Node existence can be checked by calling hasNode(u).
     */
    count degree(NodeT v) const {
        assert(hasNode(v));
        return outEdges[v].size();
    }

    /**
     * Get the number of incoming neighbors of @a v.
     *
     * @param v Node.
     * @return The number of incoming neighbors.
     * @note If the graph is not directed, the outgoing degree is returned.
     * @note The existence of the node is not checked. Calling this function with a non-existing
     * node results in a segmentation fault. Node existence can be checked by calling hasNode(u).
     */
    count degreeIn(NodeT v) const {
        assert(hasNode(v));
        return directed ? inEdges[v].size() : outEdges[v].size();
    }

    /**
     * Get the number of outgoing neighbors of @a v.
     *
     * @param v Node.
     * @return The number of outgoing neighbors.
     * @note The existence of the node is not checked. Calling this function with a non-existing
     * node results in a segmentation fault. Node existence can be checked by calling hasNode(u).
     */
    count degreeOut(NodeT v) const { return degree(v); }

    /**
     * Check whether @a v is isolated, i.e. degree is 0.
     * @param v Node.
     * @return @c true if the node is isolated (= degree is 0)
     */
    bool isIsolated(NodeT v) const {
        if (!exists[v])
            throw std::runtime_error("Error, the node does not exist!");
        return outEdges[v].empty() && (!directed || inEdges[v].empty());
    }

    /**
     * Returns the weighted degree of @a u.
     *
     * @param u Node.
     * @param countSelfLoopsTwice If set to true, self-loops will be counted twice.
     *
     * @return Weighted degree of @a u.
     */
    EdgeWeightT weightedDegree(NodeT u, bool countSelfLoopsTwice = false) const;

    /**
     * Returns the weighted in-degree of @a u.
     *
     * @param u Node.
     * @param countSelfLoopsTwice If set to true, self-loops will be counted twice.
     *
     * @return Weighted in-degree of @a v.
     */
    EdgeWeightT weightedDegreeIn(NodeT u, bool countSelfLoopsTwice = false) const;

    /* EDGE MODIFIERS */

    /**
     * Insert an edge between the nodes @a u and @a v. If the graph is
     * weighted you can optionally set a weight for this edge. The default
     * weight is 1.0. Note: Multi-edges are not supported and will NOT be
     * handled consistently by the graph data structure. It is possible to check
     * for multi-edges by enabling parameter "checkForMultiEdges". If already present,
     * the new edge is not inserted. Enabling this check increases the complexity of the function
     * to O(max(deg(u), deg(v))).
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @param ew Optional edge weight.
     * @param checkMultiEdge If true, this enables a check for a possible multi-edge.
     * @return @c true if edge has been added, false otherwise (in case checkMultiEdge is set to
     * true and the new edge would have been a multi-edge.)
     */
    bool addEdge(NodeT u, NodeT v, EdgeWeightT ew = defaultEdgeWeightT,
                 bool checkMultiEdge = false);

    /**
     * Insert an edge between the nodes @a u and @a v. Unline the addEdge function, this function
     * does not not add any information to v. If the graph is weighted you can optionally set a
     * weight for this edge. The default weight is 1.0. Note: Multi-edges are not supported and will
     * NOT be handled consistently by the graph data structure. It is possible to check
     * for multi-edges by enabling parameter "checkForMultiEdges". If already present,
     * the new edge is not inserted. Enabling this check increases the complexity of the function
     * to O(max(deg(u), deg(v))).
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @param ew Optional edge weight.
     * @param ew Optional edge weight.
     * @param index Optional edge index.
     * @param checkForMultiEdges If true, this enables a check for a possible multi-edge.
     * @return @c true if edge has been added, false otherwise (in case checkMultiEdge is set to
     * true and the new edge would have been a multi-edge.)
     */
    bool addPartialEdge(Unsafe, NodeT u, NodeT v, EdgeWeightT ew = defaultEdgeWeightT,
                        uint64_t index = 0, bool checkForMultiEdges = false);

    /**
     * Insert an in edge between the nodes @a u and @a v in a directed graph. If the graph is
     * weighted you can optionally set a weight for this edge. The default
     * weight is 1.0. Note: Multi-edges are not supported and will NOT be
     * handled consistently by the graph data structure. It is possible to check
     * for multi-edges by enabling parameter "checkForMultiEdges". If already present,
     * the new edge is not inserted. Enabling this check increases the complexity of the function
     * to O(max(deg(u), deg(v))).
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @param ew Optional edge weight.
     * @param index Optional edge index.
     * @param checkForMultiEdges If true, this enables a check for a possible multi-edge.
     * @return @c true if edge has been added, false otherwise (in case checkMultiEdge is set to
     * true and the new edge would have been a multi-edge.)
     */
    bool addPartialInEdge(Unsafe, NodeT u, NodeT v, EdgeWeightT ew = defaultEdgeWeightT,
                          uint64_t index = 0, bool checkForMultiEdges = false);

    /**
     * Insert an out edge between the nodes @a u and @a v in a directed graph. If the graph is
     * weighted you can optionally set a weight for this edge. The default
     * weight is 1.0. Note: Multi-edges are not supported and will NOT be
     * handled consistently by the graph data structure. It is possible to check
     * for multi-edges by enabling parameter "checkForMultiEdges". If already present,
     * the new edge is not inserted. Enabling this check increases the complexity of the function
     * to O(max(deg(u), deg(v))).
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @param ew Optional edge weight.
     * @param index Optional edge index.
     * @param checkForMultiEdges If true, this enables a check for a possible multi-edge.
     * @return @c true if edge has been added, false otherwise (in case checkMultiEdge is set to
     * true and the new edge would have been a multi-edge.)
     */
    bool addPartialOutEdge(Unsafe, NodeT u, NodeT v, EdgeWeightT ew = defaultEdgeWeightT,
                           uint64_t index = 0, bool checkForMultiEdges = false);

    /**
     * If set to true, the ingoing and outgoing adjacency vectors will
     * automatically be updated to maintain a sorting (if it existed before) by performing up to n-1
     * swaps. If the user plans to remove multiple edges, better set it to false and call
     * sortEdges() afterwards to avoid redundant swaps. Default = true.
     */
    void setKeepEdgesSorted(bool sorted = true) { maintainSortedEdges = sorted; }

    /**
     * If set to true, the EdgeIDs will automatically be adjusted,
     * so that no gaps in between IDs exist. If the user plans to remove multiple edges, better set
     * it to false and call indexEdges(force=true) afterwards to avoid redundant re-indexing.
     * Default = true.
     */
    void setMaintainCompactEdges(bool compact = true) { maintainCompactEdges = compact; }

    /**
     * Returns true if edges are currently being sorted when removeEdge() is called.
     */
    bool getKeepEdgesSorted() const noexcept { return maintainSortedEdges; }

    /*
     * Returns true if edges are currently being compacted when removeEdge() is called.
     */
    bool getMaintainCompactEdges() const noexcept { return maintainCompactEdges; }

    /**
     *
     * Removes the undirected edge {@a u,@a v}.
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     */
    void removeEdge(NodeT u, NodeT v);

    /**
     * Removes all the edges in the graph.
     */
    void removeAllEdges();

    /**
     * Removes edges adjacent to a node according to a specific criterion.
     *
     * @param u The node whose adjacent edges shall be removed.
     * @param condition A function that takes a node as an input and returns a
     * bool. If true the edge (u, v) is removed.
     * @param edgesIn Whether in-going or out-going edges shall be removed.
     * @return std::pair<count, count> The number of removed edges (first) and the number of removed
     * self-loops (second).
     */
    template <typename Condition>
    std::pair<count, count> removeAdjacentEdges(NodeT u, Condition condition, bool edgesIn = false);

    /**
     * Removes all self-loops in the graph.
     */
    void removeSelfLoops();

    /**
     * Removes all multi-edges in the graph.
     */
    void removeMultiEdges();

    /**
     * Changes the edges {@a s1, @a t1} into {@a s1, @a t2} and the edge {@a
     * s2,
     * @a t2} into {@a s2, @a t1}.
     *
     * If there are edge weights or edge ids, they are preserved. Note that no
     * check is performed if the swap is actually possible, i.e. does not
     * generate duplicate edges.
     *
     * @param s1 The first source
     * @param t1 The first target
     * @param s2 The second source
     * @param t2 The second target
     */
    void swapEdge(NodeT s1, NodeT t1, NodeT s2, NodeT t2);

    /**
     * Checks if undirected edge {@a u,@a v} exists in the graph.
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @return <code>true</code> if the edge exists, <code>false</code>
     * otherwise.
     */
    bool hasEdge(NodeT u, NodeT v) const noexcept;

    /* GLOBAL PROPERTIES */

    /**
     * Returns <code>true</code> if this graph supports edge weights other
     * than 1.0.
     * @return <code>true</code> if this graph supports edge weights other
     * than 1.0.
     */
    bool isWeighted() const noexcept { return weighted; }

    /**
     * Return @c true if this graph supports directed edges.
     * @return @c true if this graph supports directed edges.
     */
    bool isDirected() const noexcept { return directed; }

    /**
     * Return <code>true</code> if graph contains no nodes.
     * @return <code>true</code> if graph contains no nodes.
     */
    bool isEmpty() const noexcept { return !n; }

    /**
     * Return the number of nodes in the graph.
     * @return The number of nodes.
     */
    count numberOfNodes() const noexcept { return n; }

    /**
     * Return the number of edges in the graph.
     * @return The number of edges.
     */
    count numberOfEdges() const noexcept { return m; }

    /**
     * Return the number of loops {v,v} in the graph.
     * @return The number of loops.
     * @note This involves calculation, so store result if needed multiple
     * times.
     */
    count numberOfSelfLoops() const noexcept { return storedNumberOfSelfLoops; }

    /**
     * Get an upper bound for the node ids in the graph.
     * @return An upper bound for the node ids.
     */
    index upperNodeIdBound() const noexcept { return z; }

    /**
     * Check for invalid graph states, such as multi-edges.
     * @return False if the graph is in invalid state.
     */
    bool checkConsistency() const;

    /* DYNAMICS */

    /**
     * Trigger a time step - increments counter.
     *
     * This method is deprecated and will not be supported in future releases.
     */
    void timeStep() {
        WARN("Graph::timeStep should not be used and will be deprecated in the future.");
        ++t;
    }

    /**
     * Get time step counter.
     * @return Time step counter.
     *
     * This method is deprecated and will not be supported in future releases.
     */
    count time() {
        WARN("Graph::time should not be used and will be deprecated in the future.");
        return t;
    }

    /**
     * Return edge weight of edge {@a u,@a v}. Returns 0 if edge does not
     * exist. BEWARE: Running time is \Theta(deg(u))!
     *
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @return Edge weight of edge {@a u,@a v} or 0 if edge does not exist.
     */
    EdgeWeightT weight(NodeT u, NodeT v) const;

    /**
     * Set the weight of an edge. If the edge does not exist,
     * it will be inserted.
     *
     * @param[in]	u	endpoint of edge
     * @param[in]	v	endpoint of edge
     * @param[in]	ew	edge weight
     */
    void setWeight(NodeT u, NodeT v, EdgeWeightT ew);

    /**
     * Set the weight to the i-th neighbour of u.
     *
     * @param[in]	u	endpoint of edge
     * @param[in]	i	index of the nexight
     * @param[in]	ew	edge weight
     */
    void setWeightAtIthNeighbor(Unsafe, NodeT u, index i, EdgeWeightT ew);

    /**
     * Set the weight to the i-th incoming neighbour of u.
     *
     * @param[in]	u	endpoint of edge
     * @param[in]	i	index of the nexight
     * @param[in]	ew	edge weight
     */
    void setWeightAtIthInNeighbor(Unsafe, NodeT u, index i, EdgeWeightT ew);

    /**
     * Increase the weight of an edge. If the edge does not exist,
     * it will be inserted.
     *
     * @param[in]	u	endpoint of edge
     * @param[in]	v	endpoint of edge
     * @param[in]	ew	edge weight
     */
    void increaseWeight(NodeT u, NodeT v, EdgeWeightT ew);

    /* SUMS */

    /**
     * Returns the sum of all edge weights.
     * @return The sum of all edge weights.
     */
    EdgeWeightT totalEdgeWeight() const noexcept;

    /**
     * Return the i-th (outgoing) neighbor of @a u.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return @a i-th (outgoing) neighbor of @a u, or @c nullNodeId if no such
     * neighbor exists.
     */
    NodeT getIthNeighbor(Unsafe, NodeT u, index i) const { return outEdges[u][i]; }

    /**
     * Return the weight to the i-th (outgoing) neighbor of @a u.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return @a edge weight to the i-th (outgoing) neighbor of @a u, or @c +inf if no such
     * neighbor exists.
     */
    EdgeWeightT getIthNeighborWeight(Unsafe, NodeT u, index i) const {
        return isWeighted() ? outEdgeWeights[u][i] : defaultEdgeWeightT;
    }

    /**
     * Get an iterable range over the nodes of the graph.
     *
     * @return Iterator range over the nodes of the graph.
     */
    NodeRange nodeRange() const noexcept { return NodeRange(*this); }

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

    /**
     * Get an iterable range over the neighbors of @a.
     *
     * @param u Node.
     * @return Iterator range over the neighbors of @a u.
     */
    NeighborRange<false> neighborRange(NodeT u) const {
        assert(exists[u]);
        return NeighborRange<false>(*this, u);
    }

    /**
     * Get an iterable range over the neighbors of @a u including the edge
     * weights.
     *
     * @param u Node.
     * @return Iterator range over pairs of neighbors of @a u and corresponding
     * edge weights.
     */
    NeighborWeightRange<false> weightNeighborRange(NodeT u) const {
        assert(isWeighted());
        assert(exists[u]);
        return NeighborWeightRange<false>(*this, u);
    }

    /**
     * Get an iterable range over the in-neighbors of @a.
     *
     * @param u Node.
     * @return Iterator range over pairs of in-neighbors of @a u.
     */
    NeighborRange<true> inNeighborRange(NodeT u) const {
        assert(isDirected());
        assert(exists[u]);
        return NeighborRange<true>(*this, u);
    }

    /**
     * Get an iterable range over the in-neighbors of @a u including the
     * edge weights.
     *
     * @param u Node.
     * @return Iterator range over pairs of in-neighbors of @a u and corresponding
     * edge weights.
     */
    NeighborWeightRange<true> weightInNeighborRange(NodeT u) const {
        assert(isDirected() && isWeighted());
        assert(exists[u]);
        return NeighborWeightRange<true>(*this, u);
    }

    /**
     * Returns the index of node v in the array of outgoing edges of node u.
     *
     * @param u Node
     * @param v Node
     * @return index of node v in the array of outgoing edges of node u.
     */
    index indexOfNeighbor(NodeT u, NodeT v) const { return indexInOutEdgeArray(u, v); }

    /**
     * Return the i-th (outgoing) neighbor of @a u.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return @a i-th (outgoing) neighbor of @a u, or @c nullNodeId if no such
     * neighbor exists.
     */
    NodeT getIthNeighbor(NodeT u, index i) const {
        if (!hasNode(u) || i >= outEdges[u].size())
            return nullNodeId;
        return outEdges[u][i];
    }

    /**
     * Return the i-th (incoming) neighbor of @a u.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeIn(u))
     * @return @a i-th (incoming) neighbor of @a u, or @c nullNodeId if no such
     * neighbor exists.
     */
    NodeT getIthInNeighbor(NodeT u, index i) const {
        if (!hasNode(u) || i >= inEdges[u].size())
            return nullNodeId;
        return inEdges[u][i];
    }

    /**
     * Return the weight to the i-th (outgoing) neighbor of @a u.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return @a edge weight to the i-th (outgoing) neighbor of @a u, or @c +inf if no such
     * neighbor exists.
     */
    EdgeWeightT getIthNeighborWeight(NodeT u, index i) const {
        if (!hasNode(u) || i >= outEdges[u].size())
            return nullWeight;
        return isWeighted() ? outEdgeWeights[u][i] : defaultEdgeWeightT;
    }

    /**
     * Get i-th (outgoing) neighbor of @a u and the corresponding edge weight.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return pair: i-th (outgoing) neighbor of @a u and the corresponding
     * edge weight, or @c defaultEdgeWeightT if unweighted.
     */
    std::pair<NodeT, EdgeWeightT> getIthNeighborWithWeight(NodeT u, index i) const {
        if (!hasNode(u) || i >= outEdges[u].size())
            return {nullNodeId, EdgeWeightT{0}};
        return getIthNeighborWithWeight(unsafe, u, i);
    }

    /**
     * Get i-th (outgoing) neighbor of @a u and the corresponding edge weight.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return pair: i-th (outgoing) neighbor of @a u and the corresponding
     * edge weight, or @c defaultEdgeWeightT if unweighted.
     */
    std::pair<NodeT, EdgeWeightT> getIthNeighborWithWeight(Unsafe, NodeT u, index i) const {
        if (!isWeighted())
            return {outEdges[u][i], defaultEdgeWeightT};
        return {outEdges[u][i], outEdgeWeights[u][i]};
    }

    /**
     * Get i-th (outgoing) neighbor of @a u and the corresponding edge id.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return pair: i-th (outgoing) neighbor of @a u and the corresponding
     * edge id, or @c nullNodeId if no such neighbor exists.
     */
    std::pair<NodeT, edgeid> getIthNeighborWithId(NodeT u, index i) const {
        assert(hasEdgeIds());
        if (!hasNode(u) || i >= outEdges[u].size())
            return {nullNodeId, nullEdgeId};
        return {outEdges[u][i], outEdgeIds[u][i]};
    }

    /* NODE ITERATORS */

    /**
     * Iterate over all nodes of the graph and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L>
    void forNodes(L handle) const; // NOLINT(performance-unnecessary-value-param)

    /**
     * Iterate randomly over all nodes of the graph and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L>
    void parallelForNodes(L handle) const; // NOLINT(performance-unnecessary-value-param)

    /** Iterate over all nodes of the graph and call @a handle (lambda
     * closure) as long as @a condition remains true. This allows for breaking
     * from a node loop.
     *
     * @param condition Returning <code>false</code> breaks the loop.
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename C, typename L>
    void forNodesWhile(C condition, L handle) const; // NOLINT(performance-unnecessary-value-param)

    /**
     * Iterate randomly over all nodes of the graph and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L>
    void forNodesInRandomOrder(L handle) const; // NOLINT(performance-unnecessary-value-param)

    /**
     * Iterate in parallel over all nodes of the graph and call handler
     * (lambda closure). Using schedule(guided) to remedy load-imbalances due
     * to e.g. unequal degree distribution.
     *
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L>
    void balancedParallelForNodes(L handle) const; // NOLINT(performance-unnecessary-value-param)

    /**
     * Iterate over all undirected pairs of nodes and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameters <code>(node, node)</code>.
     */
    template <typename L>
    void forNodePairs(L handle) const; // NOLINT(performance-unnecessary-value-param)

    /**
     * Iterate over all undirected pairs of nodes in parallel and call @a
     * handle (lambda closure).
     *
     * @param handle Takes parameters <code>(node, node)</code>.
     */
    template <typename L>
    void parallelForNodePairs(L handle) const; // NOLINT(performance-unnecessary-value-param)

    /* EDGE ITERATORS */

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

    /* NEIGHBORHOOD ITERATORS */

    /**
     * Iterate over all neighbors of a node and call @a handle (lamdba
     * closure).
     *
     * @param u Node.
     * @param handle Takes parameter <code>(node)</code> or <code>(node,
     * edgeweight)</code> which is a neighbor of @a u.
     * @note For directed graphs only outgoing edges from @a u are considered.
     * A node is its own neighbor if there is a self-loop.
     *
     */
    template <typename L>
    void forNeighborsOf(NodeT u, L handle) const;

    /**
     * Iterate over all incident edges of a node and call @a handle (lamdba
     * closure).
     *
     * @param u Node.
     * @param handle Takes parameters <code>(node, node)</code>, <code>(node,
     * node, edgeweight)</code>, <code>(node, node, edgeid)</code> or
     * <code>(node, node, edgeweight, edgeid)</code> where the first node is
     * @a u and the second is a neighbor of @a u.
     * @note For undirected graphs all edges incident to @a u are also
     * outgoing edges.
     */
    template <typename L>
    void forEdgesOf(NodeT u, L handle) const;

    /**
     * Iterate over all neighbors of a node and call handler (lamdba closure).
     * For directed graphs only incoming edges from u are considered.
     */
    template <typename L>
    void forInNeighborsOf(NodeT u, L handle) const;

    /**
     * Iterate over all incoming edges of a node and call handler (lamdba
     * closure).
     * @note For undirected graphs all edges incident to u are also incoming
     * edges.
     *
     * Handle takes parameters (u, v) or (u, v, w) where w is the edge weight.
     */
    template <typename L>
    void forInEdgesOf(NodeT u, L handle) const;

    /* REDUCTION ITERATORS */

    /**
     * Iterate in parallel over all nodes and sum (reduce +) the values
     * returned by the handler
     */
    template <typename L>
    double parallelSumForNodes(L handle) const; // NOLINT(performance-unnecessary-value-param)

    /**
     * Iterate in parallel over all edges and sum (reduce +) the values
     * returned by the handler
     */
    template <typename L>
    double parallelSumForEdges(L handle) const;
};

using Graph = AdjListGraph<node, edgeweight>;

/* NODE ITERATORS */

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forNodes(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
    for (NodeT v = 0; v < z; ++v) {
        if (exists[v]) {
            handle(v);
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::parallelForNodes(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
#pragma omp parallel for
    for (omp_index v = 0; v < static_cast<omp_index>(z); ++v) {
        if (exists[v]) {
            handle(v);
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <typename C, typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forNodesWhile(
    C condition, L handle) const { // NOLINT(performance-unnecessary-value-param)
    for (NodeT v = 0; v < z; ++v) {
        if (exists[v]) {
            if (!condition()) {
                break;
            }
            handle(v);
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forNodesInRandomOrder(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
    std::vector<NodeT> randVec;
    randVec.reserve(numberOfNodes());
    forNodes([&](NodeT u) { randVec.push_back(u); });
    std::ranges::shuffle(randVec, Aux::Random::getURNG());
    for (NodeT v : randVec) {
        handle(v);
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::balancedParallelForNodes(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
// TODO: define min block size (and test it!)
#pragma omp parallel for schedule(guided)
    for (omp_index v = 0; v < static_cast<omp_index>(z); ++v) {
        if (exists[v]) {
            handle(v);
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forNodePairs(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
    for (NodeT u = 0; u < z; ++u) {
        if (exists[u]) {
            for (NodeT v = u + 1; v < z; ++v) {
                if (exists[v]) {
                    handle(u, v);
                }
            }
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::parallelForNodePairs(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
#pragma omp parallel for schedule(guided)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        if (exists[u]) {
            for (NodeT v = u + 1; v < z; ++v) {
                if (exists[v]) {
                    handle(u, v);
                }
            }
        }
    }
}

/* EDGE ITERATORS */

template <class NodeT, class EdgeWeightT>
template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline void AdjListGraph<NodeT, EdgeWeightT>::forOutEdgesOfImpl(NodeT u, L handle) const {
    for (index i = 0; i < outEdges[u].size(); ++i) {
        NodeT v = outEdges[u][i];

        if (useEdgeInIteration<graphIsDirected>(u, v)) {
            edgeLambda<L>(handle, u, v, getOutEdgeWeight<hasWeights>(u, i),
                          getOutEdgeId<graphHasEdgeIds>(u, i));
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline void AdjListGraph<NodeT, EdgeWeightT>::forInEdgesOfImpl(NodeT u, L handle) const {
    if (graphIsDirected) {
        for (index i = 0; i < inEdges[u].size(); ++i) {
            NodeT v = inEdges[u][i];

            edgeLambda<L>(handle, u, v, getInEdgeWeight<hasWeights>(u, i),
                          getInEdgeId<graphHasEdgeIds>(u, i));
        }
    } else {
        for (index i = 0; i < outEdges[u].size(); ++i) {
            NodeT v = outEdges[u][i];

            edgeLambda<L>(handle, u, v, getOutEdgeWeight<hasWeights>(u, i),
                          getOutEdgeId<graphHasEdgeIds>(u, i));
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline void AdjListGraph<NodeT, EdgeWeightT>::forEdgeImpl(L handle) const {
    for (NodeT u = 0; u < z; ++u) {
        forOutEdgesOfImpl<graphIsDirected, hasWeights, graphHasEdgeIds, L>(u, handle);
    }
}

template <class NodeT, class EdgeWeightT>
template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline void AdjListGraph<NodeT, EdgeWeightT>::parallelForEdgesImpl(L handle) const {
#pragma omp parallel for schedule(guided)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        forOutEdgesOfImpl<graphIsDirected, hasWeights, graphHasEdgeIds, L>(u, handle);
    }
}

template <class NodeT, class EdgeWeightT>
template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline double AdjListGraph<NodeT, EdgeWeightT>::parallelSumForEdgesImpl(L handle) const {
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        for (index i = 0; i < outEdges[u].size(); ++i) {
            NodeT v = outEdges[u][i];

            // undirected, do not iterate over edges twice
            // {u, v} instead of (u, v); if v == nullNodeId, u > v is not fulfilled
            if (useEdgeInIteration<graphIsDirected>(u, v)) {
                sum += edgeLambda<L>(handle, u, v, getOutEdgeWeight<hasWeights>(u, i),
                                     getOutEdgeId<graphHasEdgeIds>(u, i));
            }
        }
    }

    return sum;
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forEdges(L handle) const {
    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edgeIds
        forEdgeImpl<false, false, false, L>(handle);
        break;

    case 1: // weighted,   undirected, no edgeIds
        forEdgeImpl<false, true, false, L>(handle);
        break;

    case 2: // unweighted, directed, no edgeIds
        forEdgeImpl<true, false, false, L>(handle);
        break;

    case 3: // weighted, directed, no edgeIds
        forEdgeImpl<true, true, false, L>(handle);
        break;

    case 4: // unweighted, undirected, with edgeIds
        forEdgeImpl<false, false, true, L>(handle);
        break;

    case 5: // weighted,   undirected, with edgeIds
        forEdgeImpl<false, true, true, L>(handle);
        break;

    case 6: // unweighted, directed, with edgeIds
        forEdgeImpl<true, false, true, L>(handle);
        break;

    case 7: // weighted,   directed, with edgeIds
        forEdgeImpl<true, true, true, L>(handle);
        break;
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::parallelForEdges(L handle) const {
    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edgeIds
        parallelForEdgesImpl<false, false, false, L>(handle);
        break;

    case 1: // weighted,   undirected, no edgeIds
        parallelForEdgesImpl<false, true, false, L>(handle);
        break;

    case 2: // unweighted, directed, no edgeIds
        parallelForEdgesImpl<true, false, false, L>(handle);
        break;

    case 3: // weighted, directed, no edgeIds
        parallelForEdgesImpl<true, true, false, L>(handle);
        break;

    case 4: // unweighted, undirected, with edgeIds
        parallelForEdgesImpl<false, false, true, L>(handle);
        break;

    case 5: // weighted,   undirected, with edgeIds
        parallelForEdgesImpl<false, true, true, L>(handle);
        break;

    case 6: // unweighted, directed, with edgeIds
        parallelForEdgesImpl<true, false, true, L>(handle);
        break;

    case 7: // weighted,   directed, with edgeIds
        parallelForEdgesImpl<true, true, true, L>(handle);
        break;
    }
}

/* NEIGHBORHOOD ITERATORS */

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forNeighborsOf(NodeT u, L handle) const {
    forEdgesOf(u, handle);
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forEdgesOf(NodeT u, L handle) const {
    switch (weighted + 2 * edgesIndexed) {
    case 0: // not weighted, no edge ids
        forOutEdgesOfImpl<true, false, false, L>(u, handle);
        break;

    case 1: // weighted, no edge ids
        forOutEdgesOfImpl<true, true, false, L>(u, handle);
        break;

    case 2: // not weighted, with edge ids
        forOutEdgesOfImpl<true, false, true, L>(u, handle);
        break;

    case 3: // weighted, with edge ids
        forOutEdgesOfImpl<true, true, true, L>(u, handle);
        break;
    }
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forInNeighborsOf(NodeT u, L handle) const {
    forInEdgesOf(u, handle);
}

template <class NodeT, class EdgeWeightT>
template <typename L>
void AdjListGraph<NodeT, EdgeWeightT>::forInEdgesOf(NodeT u, L handle) const {
    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edge ids
        forInEdgesOfImpl<false, false, false, L>(u, handle);
        break;

    case 1: // weighted, undirected, no edge ids
        forInEdgesOfImpl<false, true, false, L>(u, handle);
        break;

    case 2: // unweighted, directed, no edge ids
        forInEdgesOfImpl<true, false, false, L>(u, handle);
        break;

    case 3: // weighted, directed, no edge ids
        forInEdgesOfImpl<true, true, false, L>(u, handle);
        break;

    case 4: // unweighted, undirected, with edge ids
        forInEdgesOfImpl<false, false, true, L>(u, handle);
        break;

    case 5: // weighted, undirected, with edge ids
        forInEdgesOfImpl<false, true, true, L>(u, handle);
        break;

    case 6: // unweighted, directed, with edge ids
        forInEdgesOfImpl<true, false, true, L>(u, handle);
        break;

    case 7: // weighted, directed, with edge ids
        forInEdgesOfImpl<true, true, true, L>(u, handle);
        break;
    }
}

/* REDUCTION ITERATORS */

template <class NodeT, class EdgeWeightT>
template <typename L>
double AdjListGraph<NodeT, EdgeWeightT>::parallelSumForNodes(
    L handle) const { // NOLINT(performance-unnecessary-value-param)
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (omp_index v = 0; v < static_cast<omp_index>(z); ++v) {
        if (exists[v]) {
            sum += handle(v);
        }
    }

    return sum;
}

template <class NodeT, class EdgeWeightT>
template <typename L>
double AdjListGraph<NodeT, EdgeWeightT>::parallelSumForEdges(L handle) const {
    double sum = 0.0;

    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edge ids
        sum = parallelSumForEdgesImpl<false, false, false, L>(handle);
        break;

    case 1: // weighted,   undirected, no edge ids
        sum = parallelSumForEdgesImpl<false, true, false, L>(handle);
        break;

    case 2: // unweighted, directed, no edge ids
        sum = parallelSumForEdgesImpl<true, false, false, L>(handle);
        break;

    case 3: // weighted,   directed, no edge ids
        sum = parallelSumForEdgesImpl<true, true, false, L>(handle);
        break;

    case 4: // unweighted, undirected, with edge ids
        sum = parallelSumForEdgesImpl<false, false, true, L>(handle);
        break;

    case 5: // weighted,   undirected, with edge ids
        sum = parallelSumForEdgesImpl<false, true, true, L>(handle);
        break;

    case 6: // unweighted, directed, with edge ids
        sum = parallelSumForEdgesImpl<true, false, true, L>(handle);
        break;

    case 7: // weighted,   directed, with edge ids
        sum = parallelSumForEdgesImpl<true, true, true, L>(handle);
        break;
    }

    return sum;
}

/* EDGE MODIFIERS */

template <class NodeT, class EdgeWeightT>
template <typename Condition>
std::pair<count, count>
AdjListGraph<NodeT, EdgeWeightT>::removeAdjacentEdges(NodeT u, Condition condition, bool edgesIn) {
    count removedEdges = 0;
    count removedSelfLoops = 0;

    // For directed graphs, this function is supposed to be called twice: one to remove out-edges,
    // and one to remove in-edges.
    auto& edges_ = edgesIn ? inEdges[u] : outEdges[u];
    for (index vi = 0; vi < edges_.size();) {
        if (condition(edges_[vi])) {
            auto const isSelfLoop = (edges_[vi] == u);
            removedSelfLoops += isSelfLoop;
            removedEdges += !isSelfLoop;
            edges_[vi] = edges_.back();
            edges_.pop_back();
            if (isWeighted()) {
                auto& weights_ = edgesIn ? inEdgeWeights[u] : outEdgeWeights[u];
                weights_[vi] = weights_.back();
                weights_.pop_back();
            }
            if (hasEdgeIds()) {
                auto& edgeIds_ = edgesIn ? inEdgeIds[u] : outEdgeIds[u];
                edgeIds_[vi] = edgeIds_.back();
                edgeIds_.pop_back();
            }
        } else {
            ++vi;
        }
    }

    return {removedEdges, removedSelfLoops};
}

template <class NodeT, class EdgeWeightT>
template <typename Lambda>
void AdjListGraph<NodeT, EdgeWeightT>::sortNeighbors(NodeT u, Lambda lambda) {
    if ((degreeIn(u) < 2) && (degree(u) < 2)) {
        return;
    }
    // Sort the outEdge-Attributes
    std::vector<index> outIndices(outEdges[u].size());
    std::iota(outIndices.begin(), outIndices.end(), 0);
    std::ranges::sort(outIndices,
                      [&](index a, index b) { return lambda(outEdges[u][a], outEdges[u][b]); });

    Aux::ArrayTools::applyPermutation(outEdges[u].begin(), outEdges[u].end(), outIndices.begin());

    if (weighted) {
        Aux::ArrayTools::applyPermutation(outEdgeWeights[u].begin(), outEdgeWeights[u].end(),
                                          outIndices.begin());
    }

    if (edgesIndexed) {
        Aux::ArrayTools::applyPermutation(outEdgeIds[u].begin(), outEdgeIds[u].end(),
                                          outIndices.begin());
    }

    // For directed graphs we need to sort the inEdge-Attributes separately
    if (directed) {
        std::vector<index> inIndices(inEdges[u].size());
        std::iota(inIndices.begin(), inIndices.end(), 0);

        std::ranges::sort(inIndices,
                          [&](index a, index b) { return lambda(inEdges[u][a], inEdges[u][b]); });

        Aux::ArrayTools::applyPermutation(inEdges[u].begin(), inEdges[u].end(), inIndices.begin());

        if (weighted) {
            Aux::ArrayTools::applyPermutation(inEdgeWeights[u].begin(), inEdgeWeights[u].end(),
                                              inIndices.begin());
        }

        if (edgesIndexed) {
            Aux::ArrayTools::applyPermutation(inEdgeIds[u].begin(), inEdgeIds[u].end(),
                                              inIndices.begin());
        }
    }
}

template <class NodeT, class EdgeWeightT>
template <class Lambda>
void AdjListGraph<NodeT, EdgeWeightT>::sortEdges(Lambda lambda) {

    std::vector<std::vector<index>> indicesGlobal(omp_get_max_threads());

    auto const sortAdjacencyArrays = [&](NodeT u, std::vector<NodeT>& adjList,
                                         std::vector<EdgeWeightT>& weights,
                                         std::vector<edgeid>& edgeIds) -> void {
        auto& indices = indicesGlobal[omp_get_thread_num()];
        if (adjList.size() > indices.size())
            indices.resize(adjList.size());

        const auto indicesEnd =
            indices.begin()
            + static_cast<
                std::iterator_traits<std::vector<index>::const_iterator>::difference_type>(
                adjList.size());
        std::iota(indices.begin(), indicesEnd, 0);

        if (isWeighted()) {
            if (hasEdgeIds())
                std::sort(indices.begin(), indicesEnd, [&](auto a, auto b) -> bool {
                    return lambda(WeightedEdgeWithId{u, adjList[a], weights[a], edgeIds[a]},
                                  WeightedEdgeWithId{u, adjList[b], weights[b], edgeIds[b]});
                });
            else
                std::sort(indices.begin(), indicesEnd, [&](auto a, auto b) -> bool {
                    return lambda(WeightedEdgeWithId{u, adjList[a], weights[a], 0},
                                  WeightedEdgeWithId{u, adjList[b], weights[b], 0});
                });
        } else if (hasEdgeIds())
            std::sort(indices.begin(), indicesEnd, [&](auto a, auto b) -> bool {
                return lambda(WeightedEdgeWithId{u, adjList[a], defaultEdgeWeightT, edgeIds[a]},
                              WeightedEdgeWithId{u, adjList[b], defaultEdgeWeightT, edgeIds[b]});
            });
        else
            std::sort(indices.begin(), indicesEnd, [&](auto a, auto b) -> bool {
                return lambda(WeightedEdgeWithId{u, adjList[a], defaultEdgeWeightT, 0},
                              WeightedEdgeWithId{u, adjList[b], defaultEdgeWeightT, 0});
            });

        Aux::ArrayTools::applyPermutation(adjList.begin(), adjList.end(), indices.begin());

        if (isWeighted())
            Aux::ArrayTools::applyPermutation(weights.begin(), weights.end(), indices.begin());

        if (hasEdgeIds())
            Aux::ArrayTools::applyPermutation(edgeIds.begin(), edgeIds.end(), indices.begin());
    };

    balancedParallelForNodes([&](NodeT const u) {
        if (degree(u) < 2)
            return;

        std::vector<EdgeWeightT> dummyEdgeWeights;
        std::vector<edgeid> dummyEdgeIds;
        sortAdjacencyArrays(u, outEdges[u], isWeighted() ? outEdgeWeights[u] : dummyEdgeWeights,
                            hasEdgeIds() ? outEdgeIds[u] : dummyEdgeIds);

        if (isDirected())
            sortAdjacencyArrays(u, inEdges[u], isWeighted() ? inEdgeWeights[u] : dummyEdgeWeights,
                                hasEdgeIds() ? inEdgeIds[u] : dummyEdgeIds);
    });
}

/** CONSTRUCTORS **/

template <class NodeT, class EdgeWeightT>
AdjListGraph<NodeT, EdgeWeightT>::AdjListGraph(count n, bool weighted, bool directed,
                                               bool edgesIndexed)
    : n(n), m(0), storedNumberOfSelfLoops(0), z(n), omega(0), t(0),

      weighted(weighted), // indicates whether the graph is weighted or not
      directed(directed), // indicates whether the graph is directed or not
      edgesIndexed(edgesIndexed), deletedID(nullEdgeId),
      // edges are not indexed by default

      exists(n, true),

      /* for directed graphs inEdges stores an adjacency list only considering
         incoming edges, for undirected graphs inEdges is not used*/
      inEdges(directed ? n : 0),

      /* for directed graphs outEdges stores an adjacency list only considering
      outgoing edges, for undirected graphs outEdges stores the adjacency list of
      undirected edges*/
      outEdges(n), inEdgeWeights(weighted && directed ? n : 0), outEdgeWeights(weighted ? n : 0),
      inEdgeIds(edgesIndexed && directed ? n : 0), outEdgeIds(edgesIndexed ? n : 0),
      nodeAttributeMap(this), edgeAttributeMap(this) {}

template <class NodeT, class EdgeWeightT>
AdjListGraph<NodeT, EdgeWeightT>::AdjListGraph(
    std::initializer_list<WeightedEdgeT<NodeT, EdgeWeightT>> edges)
    : AdjListGraph(0, true) {

    /* Number of nodes = highest node index + 1 */
    for (auto const& edge : edges) {
        NodeT x = std::max(edge.u, edge.v);
        while (numberOfNodes() <= x) {
            addNode();
        }
    }

    /* Now add all of the edges */
    for (auto const& edge : edges) {
        addEdge(edge.u, edge.v, edge.weight);
    }
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::preallocateUndirected(NodeT u, size_t size) {
    assert(!directed);
    assert(exists[u]);
    outEdges[u].reserve(size);
    if (weighted) {
        outEdgeWeights[u].reserve(size);
    }
    if (edgesIndexed) {
        outEdgeIds[u].reserve(size);
    }
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::preallocateDirected(NodeT u, size_t outSize, size_t inSize) {
    preallocateDirectedOutEdges(u, outSize);
    preallocateDirectedInEdges(u, inSize);
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::preallocateDirectedOutEdges(NodeT u, size_t outSize) {
    assert(directed);
    assert(exists[u]);
    outEdges[u].reserve(outSize);

    if (weighted) {
        outEdgeWeights[u].reserve(outSize);
    }
    if (edgesIndexed) {
        outEdges[u].reserve(outSize);
    }
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::preallocateDirectedInEdges(NodeT u, size_t inSize) {
    assert(directed);
    assert(exists[u]);
    inEdges[u].reserve(inSize);

    if (weighted) {
        inEdgeWeights[u].reserve(inSize);
    }
    if (edgesIndexed) {
        inEdgeIds[u].reserve(inSize);
    }
}
/** PRIVATE HELPERS **/

template <class NodeT, class EdgeWeightT>
index AdjListGraph<NodeT, EdgeWeightT>::indexInInEdgeArray(NodeT v, NodeT u) const {
    if (!directed) {
        return indexInOutEdgeArray(v, u);
    }
    for (index i = 0; i < inEdges[v].size(); ++i) {
        NodeT x = inEdges[v][i];
        if (x == u) {
            return i;
        }
    }
    return none;
}

template <class NodeT, class EdgeWeightT>
index AdjListGraph<NodeT, EdgeWeightT>::indexInOutEdgeArray(NodeT u, NodeT v) const {
    for (index i = 0; i < outEdges[u].size(); ++i) {
        NodeT x = outEdges[u][i];
        if (x == v) {
            return i;
        }
    }
    return none;
}

/** EDGE IDS **/

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::indexEdges(bool force) {
    if (edgesIndexed && !force)
        return;

    omega = 0; // reset edge ids (for re-indexing)

    outEdgeIds.clear(); // reset ids vector (for re-indexing)
    outEdgeIds.resize(outEdges.size());
    forNodes([&](NodeT u) { outEdgeIds[u].resize(outEdges[u].size(), nullEdgeId); });

    if (directed) {
        inEdgeIds.resize(inEdges.size());
        forNodes([&](NodeT u) { inEdgeIds[u].resize(inEdges[u].size(), nullEdgeId); });
    }

    // assign edge ids for edges in one direction
    forNodes([&](NodeT u) {
        for (index i = 0; i < outEdges[u].size(); ++i) {
            NodeT v = outEdges[u][i];
            if (v != nullNodeId && (directed || (u >= v))) {
                // new id
                edgeid id = omega++;
                outEdgeIds[u][i] = id;
            }
        }
    });

    // copy edge ids for the edges in the other direction. Note that
    // "indexInOutEdgeArray" is slow which is why this second loop in parallel
    // makes sense.
    if (!directed) {
        balancedParallelForNodes([&](NodeT u) {
            for (index i = 0; i < outEdges[u].size(); ++i) {
                NodeT v = outEdges[u][i];
                if (v != nullNodeId && outEdgeIds[u][i] == nullEdgeId) {
                    index j = indexInOutEdgeArray(v, u);
                    outEdgeIds[u][i] = outEdgeIds[v][j];
                }
            }
        });
    } else {
        balancedParallelForNodes([&](NodeT u) {
            for (index i = 0; i < inEdges[u].size(); ++i) {
                NodeT v = inEdges[u][i];
                if (v != nullNodeId) {
                    index j = indexInOutEdgeArray(v, u);
                    inEdgeIds[u][i] = outEdgeIds[v][j];
                }
            }
        });
    }

    edgesIndexed = true; // remember that edges have been indexed so that addEdge
                         // needs to create edge ids
}

template <class NodeT, class EdgeWeightT>
edgeid AdjListGraph<NodeT, EdgeWeightT>::edgeId(NodeT u, NodeT v) const {
    if (!edgesIndexed) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    index const i = indexInOutEdgeArray(u, v);

    if (i == none) {
        throw std::runtime_error("Edge does not exist");
    }
    return outEdgeIds[u][i];
}

/** GRAPH INFORMATION **/

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::shrinkToFit() {
    exists.shrink_to_fit();

    inEdgeWeights.shrink_to_fit();
    for (auto& w : inEdgeWeights) {
        w.shrink_to_fit();
    }

    outEdgeWeights.shrink_to_fit();
    for (auto& w : outEdgeWeights) {
        w.shrink_to_fit();
    }

    inEdges.shrink_to_fit();
    for (auto& a : inEdges) {
        a.shrink_to_fit();
    }

    outEdges.shrink_to_fit();
    for (auto& a : outEdges) {
        a.shrink_to_fit();
    }
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::compactEdges() {
    this->parallelForNodes([&](NodeT u) {
        if (degreeOut(u) == 0) {
            outEdges[u].clear();
            if (weighted)
                outEdgeWeights[u].clear();
            if (edgesIndexed)
                outEdgeIds[u].clear();
        } else {
            for (index i = 0; i < outEdges[u].size(); ++i) {
                while (i < outEdges[u].size() && outEdges[u][i] == nullNodeId) {
                    outEdges[u][i] = outEdges[u].back();
                    outEdges[u].pop_back();

                    if (weighted) {
                        outEdgeWeights[u][i] = outEdgeWeights[u].back();
                        outEdgeWeights[u].pop_back();
                    }

                    if (edgesIndexed) {
                        outEdgeIds[u][i] = outEdgeIds[u].back();
                        outEdgeIds[u].pop_back();
                    }
                }
            }
        }
        if (directed) {
            if (degreeIn(u) == 0) {
                inEdges[u].clear();
                if (weighted)
                    inEdgeWeights[u].clear();
                if (edgesIndexed)
                    inEdgeIds[u].clear();
            } else {
                for (index i = 0; i < inEdges[u].size(); ++i) {
                    while (i < inEdges[u].size() && inEdges[u][i] == nullNodeId) {
                        inEdges[u][i] = inEdges[u].back();
                        inEdges[u].pop_back();

                        if (weighted) {
                            inEdgeWeights[u][i] = inEdgeWeights[u].back();
                            inEdgeWeights[u].pop_back();
                        }

                        if (edgesIndexed) {
                            inEdgeIds[u][i] = inEdgeIds[u].back();
                            inEdgeIds[u].pop_back();
                        }
                    }
                }
            }
        }
    });
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::sortEdges() {
    std::vector<std::vector<NodeT>> targetAdjacencies(upperNodeIdBound());
    std::vector<std::vector<EdgeWeightT>> targetWeight;
    std::vector<std::vector<edgeid>> targetEdgeIds;

    if (isWeighted()) {
        targetWeight.resize(upperNodeIdBound());
        forNodes([&](NodeT u) { targetWeight[u].reserve(degree(u)); });
    }
    if (hasEdgeIds()) {
        targetEdgeIds.resize(upperNodeIdBound());
        forNodes([&](NodeT u) { targetEdgeIds[u].reserve(degree(u)); });
    }

    forNodes([&](NodeT u) { targetAdjacencies[u].reserve(degree(u)); });

    auto assignToTarget = [&](NodeT u, NodeT v, EdgeWeightT w, edgeid eid) {
        targetAdjacencies[v].push_back(u);
        if (isWeighted()) {
            targetWeight[v].push_back(w);
        }
        if (hasEdgeIds()) {
            targetEdgeIds[v].push_back(eid);
        }
    };

    forNodes([&](NodeT u) { forInEdgesOf(u, assignToTarget); });

    outEdges.swap(targetAdjacencies);
    outEdgeWeights.swap(targetWeight);
    outEdgeIds.swap(targetEdgeIds);

    if (isDirected()) {
        inEdges.swap(targetAdjacencies);
        inEdgeWeights.swap(targetWeight);
        inEdgeIds.swap(targetEdgeIds);

        forNodes([&](NodeT u) {
            targetAdjacencies[u].resize(degreeIn(u));
            targetAdjacencies[u].shrink_to_fit();
            targetAdjacencies[u].clear();
            if (isWeighted()) {
                targetWeight[u].resize(degreeIn(u));
                targetWeight[u].shrink_to_fit();
                targetWeight[u].clear();
            }
            if (hasEdgeIds()) {
                targetEdgeIds[u].resize(degreeIn(u));
                targetEdgeIds[u].shrink_to_fit();
                targetEdgeIds[u].clear();
            }
        });

        forNodes([&](NodeT u) { forEdgesOf(u, assignToTarget); });

        inEdges.swap(targetAdjacencies);
        inEdgeWeights.swap(targetWeight);
        inEdgeIds.swap(targetEdgeIds);
    }
}

template <class NodeT, class EdgeWeightT>
EdgeWeightT
AdjListGraph<NodeT, EdgeWeightT>::computeWeightedDegree(NodeT u, bool inDegree,
                                                        bool countSelfLoopsTwice) const {
    if (weighted) {
        EdgeWeightT sum{0};
        auto sumWeights = [&](NodeT v, EdgeWeightT w) {
            sum += (countSelfLoopsTwice && u == v) ? EdgeWeightT{2} * w : w;
        };
        if (inDegree) {
            forInNeighborsOf(u, sumWeights);
        } else {
            forNeighborsOf(u, sumWeights);
        }
        return sum;
    }

    count sum = inDegree ? degreeIn(u) : degreeOut(u);
    auto countSelfLoops = [&](NodeT v) { sum += (u == v); };

    if (countSelfLoopsTwice && numberOfSelfLoops()) {
        if (inDegree) {
            forInNeighborsOf(u, countSelfLoops);
        } else {
            forNeighborsOf(u, countSelfLoops);
        }
    }

    return static_cast<EdgeWeightT>(sum);
}

/** NODE MODIFIERS **/

template <class NodeT, class EdgeWeightT>
NodeT AdjListGraph<NodeT, EdgeWeightT>::addNode() {
    NodeT v = z; // node gets maximum id
    ++z;         // increment node range
    ++n;         // increment node count

    // update per node data structures
    exists.push_back(true);

    outEdges.emplace_back();
    if (weighted)
        outEdgeWeights.emplace_back();
    if (edgesIndexed)
        outEdgeIds.emplace_back();

    if (directed) {
        inEdges.emplace_back();
        if (weighted)
            inEdgeWeights.emplace_back();
        if (edgesIndexed)
            inEdgeIds.emplace_back();
    }

    return v;
}

template <class NodeT, class EdgeWeightT>
NodeT AdjListGraph<NodeT, EdgeWeightT>::addNodes(count numberOfNewNodes) {
    if (numberOfNewNodes < 10) {
        // benchmarks suggested, it's cheaper to call 10 time emplace_back than resizing.
        while (numberOfNewNodes--)
            addNode();

        return z - 1;
    }

    z += numberOfNewNodes;
    n += numberOfNewNodes;

    // update per node data structures
    exists.resize(z, true);

    outEdges.resize(z);
    if (weighted)
        outEdgeWeights.resize(z);
    if (edgesIndexed)
        outEdgeIds.resize(z);

    if (directed) {
        inEdges.resize(z);
        if (weighted)
            inEdgeWeights.resize(z);
        if (edgesIndexed)
            inEdgeIds.resize(z);
    }

    return z - 1;
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::removeNode(NodeT v) {
    assert(v < z);
    assert(exists[v]);

    // Remove all outgoing and ingoing edges
    while (!outEdges[v].empty())
        removeEdge(v, outEdges[v].front());
    if (isDirected())
        while (!inEdges[v].empty())
            removeEdge(inEdges[v].front(), v);

    // Make the attributes of this node invalid
    auto& theMap = nodeAttributeMap.attrMap;
    for (auto it = theMap.begin(); it != theMap.end(); ++it) {
        auto attributeStorageBase = it->second.get();
        attributeStorageBase->invalidate(v);
    }

    exists[v] = false;
    --n;
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::restoreNode(NodeT v) {
    assert(v < z);
    assert(!exists[v]);

    exists[v] = true;
    ++n;
}

/** NODE PROPERTIES **/

template <class NodeT, class EdgeWeightT>
EdgeWeightT AdjListGraph<NodeT, EdgeWeightT>::weightedDegree(NodeT u,
                                                             bool countSelfLoopsTwice) const {
    return computeWeightedDegree(u, false, countSelfLoopsTwice);
}

template <class NodeT, class EdgeWeightT>
EdgeWeightT AdjListGraph<NodeT, EdgeWeightT>::weightedDegreeIn(NodeT u,
                                                               bool countSelfLoopsTwice) const {
    return computeWeightedDegree(u, true, countSelfLoopsTwice);
}

/** EDGE MODIFIERS **/

template <class NodeT, class EdgeWeightT>
bool AdjListGraph<NodeT, EdgeWeightT>::addEdge(NodeT u, NodeT v, EdgeWeightT ew,
                                               bool checkForMultiEdges) {
    assert(u < z);
    assert(exists[u]);
    assert(v < z);
    assert(exists[v]);

    if (checkForMultiEdges && hasEdge(u, v)) {
        return false;
    }

    // increase number of edges
    ++m;
    outEdges[u].push_back(v);

    // if edges indexed, give new id
    if (edgesIndexed) {
        edgeid id = omega++;
        outEdgeIds[u].push_back(id);
    }

    if (directed) {
        inEdges[v].push_back(u);

        if (edgesIndexed) {
            inEdgeIds[v].push_back(omega - 1);
        }

        if (weighted) {
            inEdgeWeights[v].push_back(ew);
            outEdgeWeights[u].push_back(ew);
        }

    } else if (u == v) { // self-loop case
        if (weighted) {
            outEdgeWeights[u].push_back(ew);
        }
    } else { // undirected, no self-loop
        outEdges[v].push_back(u);

        if (weighted) {
            outEdgeWeights[u].push_back(ew);
            outEdgeWeights[v].push_back(ew);
        }

        if (edgesIndexed) {
            outEdgeIds[v].push_back(omega - 1);
        }
    }

    if (u == v) { // count self loop
        ++storedNumberOfSelfLoops;
    }

    return true;
}

template <class NodeT, class EdgeWeightT>
bool AdjListGraph<NodeT, EdgeWeightT>::addPartialEdge(Unsafe, NodeT u, NodeT v, EdgeWeightT ew,
                                                      uint64_t index, bool checkForMultiEdges) {
    assert(u < z);
    assert(exists[u]);
    assert(v < z);
    assert(exists[v]);

    if (checkForMultiEdges && (std::ranges::find(outEdges[u], v) != outEdges[u].end())) {
        return false;
    }

    outEdges[u].push_back(v);

    // if edges indexed, give new id
    if (edgesIndexed) {
        outEdgeIds[u].push_back(index);
    }
    if (weighted) {
        outEdgeWeights[u].push_back(ew);
    }

    return true;
}

template <class NodeT, class EdgeWeightT>
bool AdjListGraph<NodeT, EdgeWeightT>::addPartialOutEdge(Unsafe, NodeT u, NodeT v, EdgeWeightT ew,
                                                         uint64_t index, bool checkForMultiEdges) {
    assert(u < z);
    assert(exists[u]);
    assert(v < z);
    assert(exists[v]);

    if (checkForMultiEdges && (std::ranges::find(outEdges[u], v) != outEdges[u].end())) {
        return false;
    }

    outEdges[u].push_back(v);

    // if edges indexed, give new id
    if (edgesIndexed) {
        outEdgeIds[u].push_back(index);
    }
    if (weighted) {
        outEdgeWeights[u].push_back(ew);
    }

    return true;
}

template <class NodeT, class EdgeWeightT>
bool AdjListGraph<NodeT, EdgeWeightT>::addPartialInEdge(Unsafe, NodeT u, NodeT v, EdgeWeightT ew,
                                                        uint64_t index, bool checkForMultiEdges) {
    assert(u < z);
    assert(exists[u]);
    assert(v < z);
    assert(exists[v]);

    if (checkForMultiEdges && (std::ranges::find(inEdges[u], v) != inEdges[u].end())) {
        return false;
    }

    inEdges[u].push_back(v);

    if (edgesIndexed) {
        inEdgeIds[u].push_back(index);
    }
    if (weighted) {
        inEdgeWeights[u].push_back(ew);
    }

    return true;
}

template <class NodeT, typename T>
void erase(NodeT u, index idx, std::vector<std::vector<T>>& vec) {
    vec[u][idx] = vec[u].back();
    vec[u].pop_back();
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::removeEdge(NodeT u, NodeT v) {
    assert(u < z);
    assert(exists[u]);
    assert(v < z);
    assert(exists[v]);

    if (maintainCompactEdges && !edgesIndexed) {
        throw std::runtime_error("Edges have to be indexed if maintainCompactEdges is set to true");
    }

    index vi = indexInOutEdgeArray(u, v); // index in outEdges array
    index ui = indexInInEdgeArray(v, u);  // index in inEdges array

    if (edgesIndexed) {
        deletedID = edgeId(u, v);
    }

    if (vi == none) {
        std::stringstream strm;
        strm << "edge (" << u << "," << v << ") does not exist";
        throw std::runtime_error(strm.str());
    }

    auto const isLoop = (u == v);
    --m; // decrease number of edges
    if (isLoop)
        --storedNumberOfSelfLoops;

    // remove edge for source node
    erase(u, vi, outEdges);
    if (weighted) {
        erase(u, vi, outEdgeWeights);
    }
    if (edgesIndexed) {
        erase(u, vi, outEdgeIds);
        // Make the attributes of this edge invalid
        auto& theMap = edgeAttributeMap.attrMap;
        for (auto it = theMap.begin(); it != theMap.end(); ++it) {
            auto attributeStorageBase = it->second.get();
            attributeStorageBase->invalidate(deletedID);
        }
    }
    if (!directed && !isLoop) {
        // also remove edge for target node
        erase(v, ui, outEdges);
        if (weighted) {
            erase(v, ui, outEdgeWeights);
        }
        if (edgesIndexed) {
            erase(v, ui, outEdgeIds);
        }
    }
    if (maintainSortedEdges) {
        // initial index of deleted edge, also represents current index
        index cur = vi;

        // sort edges of source node from deleted index upwards
        while (cur + 1 < outEdges[u].size() && outEdges[u][cur] > outEdges[u][cur + 1]) {
            std::swap(outEdges[u][cur], outEdges[u][cur + 1]);
            if (edgesIndexed) {
                std::swap(outEdgeIds[u][cur], outEdgeIds[u][cur + 1]);
                // swap attributes as well
                auto& theMap = edgeAttributeMap.attrMap;
                for (auto it = theMap.begin(); it != theMap.end(); ++it) {
                    auto attributeStorageBase = it->second.get();
                    attributeStorageBase->swapData(outEdgeIds[u][cur], outEdgeIds[u][cur + 1]);
                }
            }
            ++cur;
        }

        if (!directed) {
            cur = ui;

            // sort edges of target node from deleted index upwards
            while (cur + 1 < outEdges[v].size() && outEdges[v][cur] > outEdges[v][cur + 1]) {
                std::swap(outEdges[v][cur], outEdges[v][cur + 1]);
                if (edgesIndexed) {
                    std::swap(outEdgeIds[v][cur], outEdgeIds[v][cur + 1]);
                }
                ++cur;
            }
        }
    }
    if (maintainCompactEdges) {
        // re-index edge IDs from deleted edge upwards
        balancedParallelForNodes([&](NodeT w) {
            for (index i = 0; i < outEdges[w].size(); ++i) {
                auto curID = outEdgeIds[w][i];
                if (curID > deletedID) {
                    --outEdgeIds[w][i];
                }
            }
        });
        // use erase to remove data entry at index `deletedID` and compact the data vector again
        auto& theMap = edgeAttributeMap.attrMap;
        for (auto it = theMap.begin(); it != theMap.end(); ++it) {
            auto attributeStorageBase = it->second.get();
            attributeStorageBase->erase(deletedID);
        }
    }
    if (directed) {
        assert(ui != none);

        erase(v, ui, inEdges);
        if (weighted) {
            erase(v, ui, inEdgeWeights);
        }
        if (edgesIndexed) {
            erase(v, ui, inEdgeIds);
        }
        if (maintainSortedEdges) {
            // initial index of deleted edge, also represents current index
            index cur = ui;

            // sort edges of target node from deleted index upwards
            while (cur + 1 < inEdges[v].size() && inEdges[v][cur] > inEdges[v][cur + 1]) {
                std::swap(inEdges[v][cur], inEdges[v][cur + 1]);
                if (edgesIndexed) {
                    std::swap(inEdgeIds[v][cur], inEdgeIds[v][cur + 1]);
                }
                ++cur;
            }
        }

        if (maintainCompactEdges) {
            // re-index edge ids from target node
            balancedParallelForNodes([&](NodeT w) {
                for (index i = 0; i < inEdges[w].size(); ++i) {
                    NodeT vv = inEdges[w][i];
                    if (vv != none) {
                        index j = indexInOutEdgeArray(vv, w);
                        inEdgeIds[w][i] = outEdgeIds[vv][j];
                    }
                }
            });
        }
    }
    if (maintainCompactEdges) {
        --omega; // decrease upperBound of edges
    }
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::removeAllEdges() {
    parallelForNodes([&](NodeT const u) {
        removePartialOutEdges(unsafe, u);
        if (isDirected()) {
            removePartialInEdges(unsafe, u);
        }
    });

    m = 0;
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::removeSelfLoops() {
    parallelForNodes([&](NodeT const u) {
        auto isSelfLoop = [u](const NodeT v) { return u == v; };
        removeAdjacentEdges(u, isSelfLoop);
        if (isDirected()) {
            removeAdjacentEdges(u, isSelfLoop, true);
        }
    });

    m -= storedNumberOfSelfLoops;
    storedNumberOfSelfLoops = 0;
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::removeMultiEdges() {
    count removedEdges = 0;
    count removedSelfLoops = 0;
    std::unordered_set<NodeT> nodes;

    forNodes([&](NodeT const u) {
        nodes.reserve(degree(u));
        auto isMultiedge = [&nodes](const NodeT v) { return !nodes.insert(v).second; };
        auto result = removeAdjacentEdges(u, isMultiedge);
        removedEdges += result.first;
        removedSelfLoops += result.second;
        if (isDirected()) {
            nodes.clear();
            removeAdjacentEdges(u, isMultiedge, true);
        }
        nodes.clear();
    });

    if (!isDirected()) {
        assert(!(removedEdges % 2));
        removedEdges /= 2;
    }

    m -= removedEdges + removedSelfLoops;
    storedNumberOfSelfLoops -= removedSelfLoops;
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::swapEdge(NodeT s1, NodeT t1, NodeT s2, NodeT t2) {
    index s1t1 = indexInOutEdgeArray(s1, t1);
    if (s1t1 == none)
        throw std::runtime_error("The first edge does not exist");
    index t1s1 = indexInInEdgeArray(t1, s1);

    index s2t2 = indexInOutEdgeArray(s2, t2);
    if (s2t2 == none)
        throw std::runtime_error("The second edge does not exist");
    index t2s2 = indexInInEdgeArray(t2, s2);

    std::swap(outEdges[s1][s1t1], outEdges[s2][s2t2]);

    if (directed) {
        std::swap(inEdges[t1][t1s1], inEdges[t2][t2s2]);

        if (weighted) {
            std::swap(inEdgeWeights[t1][t1s1], inEdgeWeights[t2][t2s2]);
        }

        if (edgesIndexed) {
            std::swap(inEdgeIds[t1][t1s1], inEdgeIds[t2][t2s2]);
        }
    } else {
        std::swap(outEdges[t1][t1s1], outEdges[t2][t2s2]);

        if (weighted) {
            std::swap(outEdgeWeights[t1][t1s1], outEdgeWeights[t2][t2s2]);
        }

        if (edgesIndexed) {
            std::swap(outEdgeIds[t1][t1s1], outEdgeIds[t2][t2s2]);
        }
    }
}

template <class NodeT, class EdgeWeightT>
bool AdjListGraph<NodeT, EdgeWeightT>::hasEdge(NodeT u, NodeT v) const noexcept {
    if (u >= z || v >= z) {
        return false;
    }
    if (!directed && outEdges[u].size() > outEdges[v].size()) {
        return indexInOutEdgeArray(v, u) != none;
    } else if (directed && outEdges[u].size() > inEdges[v].size()) {
        return indexInInEdgeArray(v, u) != none;
    } else {
        return indexInOutEdgeArray(u, v) != none;
    }
}

/** EDGE ATTRIBUTES **/

template <class NodeT, class EdgeWeightT>
EdgeWeightT AdjListGraph<NodeT, EdgeWeightT>::weight(NodeT u, NodeT v) const {
    index vi = indexInOutEdgeArray(u, v);
    if (vi == none) {
        return nullWeight;
    } else {
        return weighted ? outEdgeWeights[u][vi] : defaultEdgeWeightT;
    }
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::setWeight(NodeT u, NodeT v, EdgeWeightT ew) {
    if (!weighted) {
        throw std::runtime_error("Cannot set edge weight in unweighted graph.");
    }

    index vi = indexInOutEdgeArray(u, v);
    if (vi == none) {
        // edge does not exist, create it, but warn user
        TRACE("Setting edge weight of a nonexisting edge will create the edge.");
        addEdge(u, v, ew);
        return;
    }

    outEdgeWeights[u][vi] = ew;
    if (directed) {
        index ui = indexInInEdgeArray(v, u);
        inEdgeWeights[v][ui] = ew;
    } else if (u != v) {
        index ui = indexInInEdgeArray(v, u);
        outEdgeWeights[v][ui] = ew;
    }
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::increaseWeight(NodeT u, NodeT v, EdgeWeightT ew) {
    if (!weighted) {
        throw std::runtime_error("Cannot increase edge weight in unweighted graph.");
    }

    index vi = indexInOutEdgeArray(u, v);
    if (vi == none) {
        // edge does not exits, create it, but warn user
        addEdge(u, v, ew);
        return;
    }

    outEdgeWeights[u][vi] += ew;
    if (directed) {
        index ui = indexInInEdgeArray(v, u);
        inEdgeWeights[v][ui] += ew;
    } else if (u != v) {
        index ui = indexInInEdgeArray(v, u);
        outEdgeWeights[v][ui] += ew;
    }
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::setWeightAtIthNeighbor(Unsafe, NodeT u, index i,
                                                              EdgeWeightT ew) {
    outEdgeWeights[u][i] = ew;
}

template <class NodeT, class EdgeWeightT>
void AdjListGraph<NodeT, EdgeWeightT>::setWeightAtIthInNeighbor(Unsafe, NodeT u, index i,
                                                                EdgeWeightT ew) {
    inEdgeWeights[u][i] = ew;
}

/** SUMS **/

template <class NodeT, class EdgeWeightT>
EdgeWeightT AdjListGraph<NodeT, EdgeWeightT>::totalEdgeWeight() const noexcept {
    if (weighted)
        return parallelSumForEdges([](NodeT, NodeT, EdgeWeightT ew) { return ew; });
    return numberOfEdges() * defaultEdgeWeightT;
}

template <class NodeT, class EdgeWeightT>
bool AdjListGraph<NodeT, EdgeWeightT>::checkConsistency() const {
    // check for multi-edges
    std::vector<NodeT> lastSeen(z, nullNodeId);
    bool noMultiEdges = true;
    auto noMultiEdgesDetected = [&noMultiEdges]() { return noMultiEdges; };
    forNodesWhile(noMultiEdgesDetected, [&](NodeT v) {
        forNeighborsOf(v, [&](NodeT u) {
            if (lastSeen[u] == v) {
                noMultiEdges = false;
                DEBUG("Multiedge found between ", u, " and ", v, "!");
            }
            lastSeen[u] = v;
        });
    });

    bool correctNodeUpperbound = (z == outEdges.size()) && ((directed ? z : 0) == inEdges.size())
                                 && ((weighted ? z : 0) == outEdgeWeights.size())
                                 && ((weighted && directed ? z : 0) == inEdgeWeights.size())
                                 && ((edgesIndexed ? z : 0) == outEdgeIds.size())
                                 && ((edgesIndexed && directed ? z : 0) == inEdgeIds.size());

    if (!correctNodeUpperbound)
        DEBUG("Saved node upper bound doesn't actually match the actual node upper bound!");

    count NumberOfOutEdges = 0;
    count NumberOfOutEdgeWeights = 0;
    count NumberOfOutEdgeIds = 0;
    for (index i = 0; i < outEdges.size(); ++i) {
        NumberOfOutEdges += outEdges[i].size();
    }
    if (weighted)
        for (index i = 0; i < outEdgeWeights.size(); ++i) {
            NumberOfOutEdgeWeights += outEdgeWeights[i].size();
        }
    if (edgesIndexed)
        for (index i = 0; i < outEdgeIds.size(); ++i) {
            NumberOfOutEdgeIds += outEdgeIds[i].size();
        }

    count NumberOfInEdges = 0;
    count NumberOfInEdgeWeights = 0;
    count NumberOfInEdgeIds = 0;
    if (directed) {
        for (index i = 0; i < inEdges.size(); ++i) {
            NumberOfInEdges += inEdges[i].size();
        }
        if (weighted)
            for (index i = 0; i < inEdgeWeights.size(); ++i) {
                NumberOfInEdgeWeights += inEdgeWeights[i].size();
            }
        if (edgesIndexed)
            for (index i = 0; i < inEdgeIds.size(); ++i) {
                NumberOfInEdgeIds += inEdgeIds[i].size();
            }
    }

    if (!directed) {
        NumberOfOutEdges = (NumberOfOutEdges + storedNumberOfSelfLoops) / 2;
        if (weighted)
            NumberOfOutEdgeWeights = (NumberOfOutEdgeWeights + storedNumberOfSelfLoops) / 2;
        if (edgesIndexed)
            NumberOfOutEdgeIds = (NumberOfOutEdgeIds + storedNumberOfSelfLoops) / 2;
    }

    bool correctNumberOfEdges = (m == NumberOfOutEdges) && ((directed ? m : 0) == NumberOfInEdges)
                                && ((weighted ? m : 0) == NumberOfOutEdgeWeights)
                                && ((weighted && directed ? m : 0) == NumberOfInEdgeWeights)
                                && ((edgesIndexed ? m : 0) == NumberOfOutEdgeIds)
                                && ((edgesIndexed && directed ? m : 0) == NumberOfInEdgeIds);

    if (!correctNumberOfEdges)
        DEBUG("Saved number of edges is incorrect!");

    return noMultiEdges && correctNodeUpperbound && correctNumberOfEdges;
}

} /* namespace NetworKit */

#endif // NETWORKIT_GRAPH_GRAPH_HPP_
