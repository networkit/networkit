#ifndef NETWORKIT_GRAPH_DHB_GRAPH_HPP_
#define NETWORKIT_GRAPH_DHB_GRAPH_HPP_

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/ArrayTools.hpp>
#include <networkit/auxiliary/FunctionTraits.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/Attributes.hpp>
#include <networkit/graph/Edge.hpp>
#include <networkit/graph/EdgeIterators.hpp>
#include <networkit/graph/NodeIterators.hpp>

#include <omp.h>
#include <dhb/dynamic_hashed_blocks.h>
#include <tlx/define/deprecated.hpp>

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <queue>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <thread>
#include <typeindex>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace NetworKit {

// forward declaration to randomization/CurveballImpl.hpp
namespace CurveballDetails {
class CurveballMaterialization;
}

/**
 * @ingroup graph
 * A graph (with optional weights) and parallel iterator methods.
 *
 * TODO:
 * - [ ] Review the passing by value semantics for all lambdas passed
 *       to a method such as f(L handle). Instead, we might want to pass
 *       all lambdas by rvalue-reference as in f(L&& handle).
 */
class DHBGraph final {

    // DHBGraph properties
    //!< true if the graph is weighted, false otherwise
    bool weighted;
    //!< true if the graph is directed, false otherwise
    bool directed;

    // Edge properties
    //!< current number of edges
    count m;
    //!< current number of self loops, edges which have the same origin and
    //!< target
    count storedNumberOfSelfLoops;
    //!< current upper bound of edge ids, will be the id of the next edge
    edgeid omega;
    //!< true if edge ids have been assigned
    bool edgesIndexed;
    struct EdgeData {
        edgeweight weight;
        edgeid id;
    };

private:
    bool addEdge(node u, node v, edgeweight ew, edgeid id);

    AttributeMap<PerNode, DHBGraph> nodeAttributeMap;
    AttributeMap<PerEdge, DHBGraph> edgeAttributeMap;
    dhb::Matrix<EdgeData> m_dhb_graph;

public:
    auto &nodeAttributes() noexcept { return nodeAttributeMap; }
    const auto &nodeAttributes() const noexcept { return nodeAttributeMap; }
    auto &edgeAttributes() noexcept { return edgeAttributeMap; }
    const auto &edgeAttributes() const noexcept { return edgeAttributeMap; }

    // wrap up some typed attributes for the cython interface:
    //

    auto attachNodeIntAttribute(const std::string &name) {
        nodeAttributes().theGraph = this;
        return nodeAttributes().attach<int>(name);
    }

    auto attachEdgeIntAttribute(const std::string &name) {
        edgeAttributes().theGraph = this;
        return edgeAttributes().attach<int>(name);
    }

    auto attachNodeDoubleAttribute(const std::string &name) {
        nodeAttributes().theGraph = this;
        return nodeAttributes().attach<double>(name);
    }

    auto attachEdgeDoubleAttribute(const std::string &name) {
        edgeAttributes().theGraph = this;
        return edgeAttributes().attach<double>(name);
    }

    auto attachNodeStringAttribute(const std::string &name) {
        nodeAttributes().theGraph = this;
        return nodeAttributes().attach<std::string>(name);
    }

    auto attachEdgeStringAttribute(const std::string &name) {
        edgeAttributes().theGraph = this;
        return edgeAttributes().attach<std::string>(name);
    }

    auto getNodeIntAttribute(const std::string &name) {
        nodeAttributes().theGraph = this;
        return nodeAttributes().get<int>(name);
    }

    auto getEdgeIntAttribute(const std::string &name) {
        edgeAttributes().theGraph = this;
        return edgeAttributes().get<int>(name);
    }

    auto getNodeDoubleAttribute(const std::string &name) {
        nodeAttributes().theGraph = this;
        return nodeAttributes().get<double>(name);
    }

    auto getEdgeDoubleAttribute(const std::string &name) {
        edgeAttributes().theGraph = this;
        return edgeAttributes().get<double>(name);
    }

    auto getNodeStringAttribute(const std::string &name) {
        nodeAttributes().theGraph = this;
        return nodeAttributes().get<std::string>(name);
    }

    auto getEdgeStringAttribute(const std::string &name) {
        edgeAttributes().theGraph = this;
        return edgeAttributes().get<std::string>(name);
    }

    void detachNodeAttribute(std::string const &name) {
        nodeAttributes().theGraph = this;
        nodeAttributes().detach(name);
    }

    void detachEdgeAttribute(std::string const &name) {
        edgeAttributes().theGraph = this;
        edgeAttributes().detach(name);
    }

    using NodeIntAttribute = Attribute<PerNode, DHBGraph, int, false>;
    using NodeDoubleAttribute = Attribute<PerNode, DHBGraph, double, false>;
    using NodeStringAttribute = Attribute<PerNode, DHBGraph, std::string, false>;

    using EdgeIntAttribute = Attribute<PerEdge, DHBGraph, int, false>;
    using EdgeDoubleAttribute = Attribute<PerEdge, DHBGraph, double, false>;
    using EdgeStringAttribute = Attribute<PerEdge, DHBGraph, std::string, false>;

private:
    /**
     * Computes the weighted in/out degree of node @a u. If graph is directed and @a inDegree is
     * true, the function runs in O(n) since it needs to iterate over all edges.
     *
     * @param u Node.
     * @param inDegree whether to compute the in degree or the out degree.
     * @param countSelfLoopsTwice If set to true, self-loops will be counted twice.
     *
     * @return Weighted in/out degree of node @a u.
     */
    edgeweight computeWeightedDegree(node u, bool inDegree = false,
                                     bool countSelfLoopsTwice = false) const;

    /**
     * This is a private helper function of dhb.
     * @param u Node.
     * @param i The index of the incoming neighbor to check.
     *
     * @return bool `true` if the node exists and the index is within the valid range of incoming
     * neighbors.
     */
    bool ithNeighborExists(node u, index i) const {
        assert(u < upperNodeIdBound());
        return i < m_dhb_graph.neighbors(u).degree();
    }

    /**
     * Removes out-going edges from node @u.
     *
     * @param node u Node.
     */
    void removePartialOutEdges(Unsafe, node u) {
        assert(u < upperNodeIdBound());
        m_dhb_graph.neighbors(u).clear();
    }

    /**
     * Removes in-going edges to node @v.
     *
     * @param node v Node.
     */
    void removePartialInEdges(Unsafe, node v) {
        assert(v < upperNodeIdBound());

        auto remove_edge = [&](dhb::Vertex const u) {
            if (m_dhb_graph.neighbors(u).exists(v)) {
                removeEdge(u, v);
            }
        };

        m_dhb_graph.for_nodes(std::move(remove_edge));
    }

    /**
     * @brief Returns if the edge (u, v) shall be used in the iteration of all
     * edgesIndexed
     *
     * @param u The source node of the edge
     * @param v The target node of the edge
     * @return If the node shall be used, i.e. if v is not none and in the
     * undirected case if u >= v
     */
    template <bool graphIsDirected>
    inline bool useEdgeInIteration(node u, node v) const;

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
    void forOutEdgesOfImpl(node u, L handle) const;

    /**
     * @brief Get edge weight and edge ID for a DHB edge
     */
    template <bool hasWeights, bool graphHasEdgeIds>
    std::tuple<edgeweight, edgeid> getDHBEdgeData(node u, node v) const;

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
    void forOutEdgesOfImplParallel(node u, L handle) const;

    /**
     * @brief Implementation of the for loop for incoming edges of u
     *
     * For undirected graphs, this is the same as forOutEdgesOfImpl but u and v
     * are changed in the handle
     * For directed graph, the function runs in O(n) since it needs to iterate over all edges to
     * find InEdges.
     * TODO: Find a way to parallelize both for loop with `collapse`
     * @param u The node
     * @param handle The handle that shall be executed for each edge
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
    inline void forInEdgesOfImpl(node u, L handle) const;

    /**
     * @brief Parallel version of `forInEdgesOfImpl`
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
    inline void forInEdgesOfImplParallel(node u, L handle) const;

    /**
     * @brief Implementation of the for loop for all edges, @see forEdges
     *
     * @param handle The handle that shall be executed for all edges
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
    void forEdgeImpl(L handle) const;

    /**
     * @brief Parallel implementation of the for loop for all edges, @see
     * parallelForEdges
     *
     * @param handle The handle that shall be executed for all edges
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
    void forEdgesImplParallel(L handle) const;

    /**
     * @brief Summation variant of the parallel for loop for all edges, @see
     * parallelSumForEdges
     *
     * @param handle The handle that shall be executed for all edges
     * @return void
     */
    template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
    inline double parallelSumForEdgesImpl(L handle) const;

    /**
     * Iterate in parallel over all nodes of the graph and call handler
     * (lambda closure). Using schedule(guided) to remedy load-imbalances due
     * to e.g. unequal degree distribution.
     *
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L>
    void balancedParallelForNodes(L handle) const;

    /*
     * In the following definition, Aux::FunctionTraits is used in order to only
     * execute lambda functions with the appropriate parameters. The
     * decltype-return type is used for determining the return type of the
     * lambda (needed for summation) but also determines if the lambda accepts
     * the correct number of parameters. Otherwise the return type declaration
     * fails and the function is excluded from overload resolution. Then there
     * are multiple possible lambdas with three (third parameter id or weight)
     * and two (second parameter can be second node id or edge weight for
     * neighbor iterators). This is checked using Aux::FunctionTraits and
     * std::enable_if. std::enable_if only defines the type member when the
     * given bool is true, this bool comes from std::is_same which compares two
     * types. The function traits give either the parameter type or if it is out
     * of bounds they define type as void.
     */

    /**
     * Triggers a static assert error when no other method is chosen. Because of
     * the use of "..." as arguments, the priority of this method is lower than
     * the priority of the other methods. This method avoids ugly and unreadable
     * template substitution error messages from the other declarations.
     */
    template <class F, void * = (void *)0>
    typename Aux::FunctionTraits<F>::result_type edgeLambda(F &, ...) const {
        // the strange condition is used in order to delay the evaluation of the
        // static assert to the moment when this function is actually used
        static_assert(!std::is_same<F, F>::value,
                      "Your lambda does not support the required parameters or the "
                      "parameters have the wrong type.");
        return std::declval<typename Aux::FunctionTraits<F>::result_type>(); // use the correct
                                                                             // return type (this
                                                                             // won't compile)
    }

    /**
     * Calls the given function f if its fourth argument is of the type edgeid
     * and third of type edgeweight Note that the decltype check is not enough
     * as edgeweight can be casted to node and we want to assure that .
     */
    template <class F,
              typename std::enable_if<
                  (Aux::FunctionTraits<F>::arity >= 3)
                  && std::is_same<edgeweight,
                                  typename Aux::FunctionTraits<F>::template arg<2>::type>::value
                  && std::is_same<edgeid, typename Aux::FunctionTraits<F>::template arg<3>::type>::
                      value>::type * = (void *)0>
    auto edgeLambda(F &f, node u, node v, edgeweight ew, edgeid id) const
        -> decltype(f(u, v, ew, id)) {
        return f(u, v, ew, id);
    }

    /**
     * Calls the given function f if its third argument is of the type edgeid,
     * discards the edge weight Note that the decltype check is not enough as
     * edgeweight can be casted to node.
     */
    template <
        class F,
        typename std::enable_if<
            (Aux::FunctionTraits<F>::arity >= 2)
            && std::is_same<edgeid, typename Aux::FunctionTraits<F>::template arg<2>::type>::value
            && std::is_same<node, typename Aux::FunctionTraits<F>::template arg<1>::type>::
                value /* prevent f(v, weight, eid)
                       */
            >::type * = (void *)0>
    auto edgeLambda(F &f, node u, node v, edgeweight, edgeid id) const -> decltype(f(u, v, id)) {
        return f(u, v, id);
    }

    /**
     * Calls the given function f if its third argument is of type edgeweight,
     * discards the edge id Note that the decltype check is not enough as node
     * can be casted to edgeweight.
     */
    template <class F,
              typename std::enable_if<
                  (Aux::FunctionTraits<F>::arity >= 2)
                  && std::is_same<edgeweight, typename Aux::FunctionTraits<F>::template arg<
                                                  2>::type>::value>::type * = (void *)0>
    auto edgeLambda(F &f, node u, node v, edgeweight ew, edgeid /*id*/) const
        -> decltype(f(u, v, ew)) {
        return f(u, v, ew);
    }

    /**
     * Calls the given function f if it has only two arguments and the second
     * argument is of type node, discards edge weight and id Note that the
     * decltype check is not enough as edgeweight can be casted to node.
     */
    template <class F, typename std::enable_if<
                           (Aux::FunctionTraits<F>::arity >= 1)
                           && std::is_same<node, typename Aux::FunctionTraits<F>::template arg<
                                                     1>::type>::value>::type * = (void *)0>
    auto edgeLambda(F &f, node u, node v, edgeweight /*ew*/, edgeid /*id*/) const
        -> decltype(f(u, v)) {
        return f(u, v);
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
    auto edgeLambda(F &f, node, node v, edgeweight ew, edgeid /*id*/) const -> decltype(f(v, ew)) {
        return f(v, ew);
    }

    /**
     * Calls the given function f if it has only one argument, discards the
     * first node id, the edge weight and the edge id
     */
    template <class F, void * = (void *)0>
    auto edgeLambda(F &f, node, node v, edgeweight, edgeid) const -> decltype(f(v)) {
        return f(v);
    }

    /**
     * Calls the given BFS handle with distance parameter
     */
    template <class F>
    auto callBFSHandle(F &f, node u, count dist) const -> decltype(f(u, dist)) {
        return f(u, dist);
    }

    /**
     * Calls the given BFS handle without distance parameter
     */
    template <class F>
    auto callBFSHandle(F &f, node u, count) const -> decltype(f(u)) {
        return f(u);
    }

public:
    // For support of API: NetworKit::DHBGraph::NodeIterator
    using NodeIterator = NodeIteratorBase<DHBGraph>;
    // For support of API: NetworKit::DHBGraph::NodeRange
    using NodeRange = NodeRangeBase<DHBGraph>;

    // For support of API: NetworKit::DHBGraph:EdgeIterator
    using EdgeIterator = EdgeTypeIterator<DHBGraph, Edge>;
    // For support of API: NetworKit::DHBGraph:EdgeWeightIterator
    using EdgeWeightIterator = EdgeTypeIterator<DHBGraph, WeightedEdge>;
    // For support of API: NetworKit::DHBGraph:EdgeRange
    using EdgeRange = EdgeTypeRange<DHBGraph, Edge>;
    // For support of API: NetworKit::DHBGraph:EdgeWeightRange
    using EdgeWeightRange = EdgeTypeRange<DHBGraph, WeightedEdge>;

    /**
     * Class to iterate over the in/out neighbors of a node.
     */
    class NeighborIterator {
        using IteratorType = dhb::BlockState<EdgeData>::const_iterator;
        IteratorType it;

    public:
        // The value type of the neighbors (i.e. nodes). Returned by
        // operator*().
        using value_type = node;

        // Reference to the value_type, required by STL.
        using reference = value_type &;

        // Pointer to the value_type, required by STL.
        using pointer = value_type *;

        // STL iterator category.
        using iterator_category = std::forward_iterator_tag;

        // Signed integer type of the result of subtracting two pointers,
        // required by STL.
        using difference_type = ptrdiff_t;

        // Own type.
        using self = NeighborIterator;

        NeighborIterator(IteratorType it) : it(it) {}

        /**
         * @brief WARNING: This contructor is required for Python and should not be used as the
         * iterator is not initialized.
         */
        NeighborIterator() {}

        NeighborIterator &operator++() {
            ++it;
            return *this;
        }

        NeighborIterator operator++(int) {
            auto const tmp = *this;
            ++(*this);
            return tmp;
        }

        NeighborIterator operator--() {
            --it;
            return *this;
        }

        NeighborIterator operator--(int) {
            auto const tmp = *this;
            --(*this);
            return tmp;
        }

        bool operator==(NeighborIterator const &rhs) const { return it == rhs.it; }

        bool operator!=(NeighborIterator const &rhs) const { return !(it == rhs.it); }

        node operator*() const { return it->vertex(); }
    };

    /**
     * Class to iterate over the in/out neighbors of a node including the edge
     * weights. Values are std::pair<node, edgeweight>.
     */
    class NeighborWeightIterator {
        using IteratorType = dhb::BlockState<EdgeData>::const_iterator;
        IteratorType it;

    public:
        // The value type of the neighbors (i.e. nodes). Returned by
        // operator*().
        using value_type = std::pair<node, edgeweight>;

        // Reference to the value_type, required by STL.
        using reference = value_type &;

        // Pointer to the value_type, required by STL.
        using pointer = value_type *;

        // STL iterator category.
        using iterator_category = std::forward_iterator_tag;

        // Signed integer type of the result of subtracting two pointers,
        // required by STL.
        using difference_type = ptrdiff_t;

        // Own type.
        using self = NeighborWeightIterator;

        NeighborWeightIterator(IteratorType it) : it(it) {}

        /**
         * @brief WARNING: This contructor is required for Python and should not be used as the
         * iterator is not initialized.
         */
        NeighborWeightIterator() {}

        NeighborWeightIterator &operator++() {
            ++it;
            return *this;
        }

        NeighborWeightIterator operator++(int) {
            auto const tmp = *this;
            ++(*this);
            return tmp;
        }

        NeighborWeightIterator operator--() {
            --it;
            return *this;
        }

        NeighborWeightIterator operator--(int) {
            auto const tmp = *this;
            --(*this);
            return tmp;
        }

        bool operator==(NeighborWeightIterator const &rhs) const { return it == rhs.it; }

        bool operator!=(NeighborWeightIterator const &rhs) const { return !(it == rhs.it); }

        std::pair<node, edgeweight> operator*() const {
            return std::make_pair(it->vertex(), it->data().weight);
        }
    };

    /**
     * Wrapper class to iterate over a range of the neighbors of a node within
     * a for loop.
     */
    template <bool InEdges = false>
    class NeighborRange {
        DHBGraph const *G;
        node u;

    public:
        NeighborRange(DHBGraph const &G, node u) : G(&G), u(u) { assert(G.hasNode(u)); };

        NeighborRange() : G(nullptr){};

        NeighborIterator begin() const {
            assert(G);
            // in DHB we don't support iterator over InEdges = true;
            static_assert(!InEdges);
            return NeighborIterator(G->m_dhb_graph.neighbors(u).begin());
        }

        NeighborIterator end() const {
            // in DHB we don't support iterator over InEdges = true;
            static_assert(!InEdges);
            return NeighborIterator(G->m_dhb_graph.neighbors(u).end());
        }
    };

    using OutNeighborRange = NeighborRange<false>;
    /**
     * Wrapper class to iterate over a range of the neighbors of a node
     * including the edge weights within a for loop.
     * Values are std::pair<node, edgeweight>.
     */
    template <bool InEdges = false>
    class NeighborWeightRange {

        const DHBGraph *G;
        node u;

    public:
        NeighborWeightRange(const DHBGraph &G, node u) : G(&G), u(u) { assert(G.hasNode(u)); };

        NeighborWeightRange() : G(nullptr){};

        NeighborWeightIterator begin() const {
            assert(G);
            // in DHB we don't support iterator over InEdges = true;
            static_assert(!InEdges);
            return NeighborWeightIterator(G->m_dhb_graph.neighbors(u).begin());
        }

        NeighborWeightIterator end() const {
            assert(G);
            // in DHB we don't support iterator over InEdges = true;
            static_assert(!InEdges);
            return NeighborWeightIterator(G->m_dhb_graph.neighbors(u).end());
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
    DHBGraph(count n = 0, bool weighted = false, bool directed = false, bool edgesIndexed = false);

    template <class EdgeMerger = std::plus<edgeweight>>
    DHBGraph(DHBGraph const &G, bool weighted, bool directed, bool edgesIndexed = false)
        : m(G.m), storedNumberOfSelfLoops(G.storedNumberOfSelfLoops),
          omega(edgesIndexed ? G.omega : 0), weighted(weighted), directed(directed),
          edgesIndexed(edgesIndexed), // edges are not indexed by default

          // empty node attribute map as last member for this graph
          nodeAttributeMap(this), edgeAttributeMap(this), m_dhb_graph(G.m_dhb_graph) {
        // G - The original
        // this - The copy

        bool const copy_weighted_graph_to_unweighted = G.isWeighted() && !weighted;
        bool const copy_directed_graph_to_undirected = G.isDirected() && !directed;

        if (copy_weighted_graph_to_unweighted) {
#pragma omp parallel for schedule(guided)
            for (omp_index u = 0; u < static_cast<omp_index>(m_dhb_graph.vertices_count()); ++u) {
                dhb::Matrix<EdgeData>::NeighborView n = m_dhb_graph.neighbors(u);
                for (auto v = n.begin(); v != n.end(); ++v) {
                    v->data().weight = defaultEdgeWeight;
                }
            }
        }

        if (copy_directed_graph_to_undirected) { // insert bidirectional edges if they don't exist
            for (dhb::Vertex u = 0u; u < m_dhb_graph.vertices_count(); ++u) {
                auto n = m_dhb_graph.neighbors(u);
                for (auto v = n.begin(); v != n.end(); ++v) {
                    if (!hasEdge(v->vertex(), u)) {
                        addEdge(v->vertex(), u, v->data().weight);
                    }
                }
            }
        }
    }

    /**
     * Generate a weighted graph from a list of edges. (Useful for small
     * graphs in unit tests that you do not want to read from a file.)
     *
     * @param[in] edges list of weighted edges
     */
    DHBGraph(std::initializer_list<WeightedEdge> edges);

    /**
     * Create a graph as copy of @a other.
     * @param other The graph to copy.
     */
    DHBGraph(const DHBGraph &other)
        : weighted(other.weighted), directed(other.directed), m(other.m),
          storedNumberOfSelfLoops(other.storedNumberOfSelfLoops), omega(other.omega),
          edgesIndexed(other.edgesIndexed), nodeAttributeMap(other.nodeAttributeMap, this),
          edgeAttributeMap(other.edgeAttributeMap, this), m_dhb_graph(other.m_dhb_graph){};

    /** Default move constructor */
    DHBGraph(DHBGraph &&other) noexcept
        : weighted(other.weighted), directed(other.directed), m(other.m),
          storedNumberOfSelfLoops(other.storedNumberOfSelfLoops), omega(other.omega),
          edgesIndexed(other.edgesIndexed), nodeAttributeMap(std::move(other.nodeAttributeMap)),
          edgeAttributeMap(std::move(other.edgeAttributeMap)), m_dhb_graph(other.m_dhb_graph) {
        // attributes: set graph pointer to this new graph
        nodeAttributeMap.theGraph = this;
        edgeAttributeMap.theGraph = this;
    };

    /** Default destructor */
    ~DHBGraph() = default;

    /** Default move assignment operator */
    DHBGraph &operator=(DHBGraph &&other) noexcept {
        std::swap(weighted, other.weighted);
        std::swap(directed, other.directed);
        std::swap(m, other.m);
        std::swap(storedNumberOfSelfLoops, other.storedNumberOfSelfLoops);
        std::swap(omega, other.omega);
        std::swap(edgesIndexed, other.edgesIndexed);

        // attributes: set graph pointer to this new graph
        std::swap(nodeAttributeMap, other.nodeAttributeMap);
        std::swap(edgeAttributeMap, other.edgeAttributeMap);
        nodeAttributeMap.theGraph = this;
        edgeAttributeMap.theGraph = this;

        std::swap(m_dhb_graph, other.m_dhb_graph);

        return *this;
    };

    /** Default copy assignment operator */
    DHBGraph &operator=(const DHBGraph &other) {
        weighted = other.weighted;
        directed = other.directed;
        m = other.m;
        storedNumberOfSelfLoops = other.storedNumberOfSelfLoops;
        omega = other.omega;
        edgesIndexed = other.edgesIndexed;

        // call special constructors to copy attribute maps
        nodeAttributeMap = AttributeMap(other.nodeAttributeMap, this);
        edgeAttributeMap = AttributeMap(other.edgeAttributeMap, this);

        m_dhb_graph = other.m_dhb_graph;

        return *this;
    };

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
    edgeid edgeId(node u, node v) const;

    /**
     * Get the Edge (u,v) of the given id. (inverse to edgeId)
     * @note Time complexity of this function is O(n).
     */
    std::pair<node, node> edgeById(index id) const {
        std::pair<node, node> result{none, none};
        bool found = false;

        forNodesWhile([&] { return !found; },
                      [&](node u) {
                          forNeighborsOf(u, [&](node v) {
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
     * Sort edges according to accending neighbors for all vertices.
     */
    void sortEdges() {
        return sortEdges([](node a_vertex, edgeweight, edgeid, node b_vertex, edgeweight, edgeid) {
            return a_vertex < b_vertex;
        });
    }
    /**
     * Default graph class: Sorts the adjacency arrays by a custom criterion.
     * DHB: Sorts the neighbors by a custom criterion.
     * @param lambda Lambda function used to sort the edges.
     * Default graph class: It takes two WeightedEdge e1 and e2 as input parameters, returns true if
     * e1 < e2, false otherwise.
     * DHB: It takes node a_vertex, edgeweight a_weight, edgeid a_id, node b_vertex, edgeweight
     * b_weight, edgeid b_id,  as input parameters, and compares them according to lambda.
     * TODO: Maybe we need to think of a nicer way for the interface in lambda.
     */
    template <typename Lambda>
    void sortEdges(Lambda &&lambda) {
        for (node v = 0; v < m_dhb_graph.vertices_count(); ++v) {
            m_dhb_graph.sort(v, [l = std::forward<Lambda>(lambda)](
                                    dhb::BlockState<EdgeData>::Entry &a,
                                    dhb::BlockState<EdgeData>::Entry &b) {
                return l(a.vertex, a.data.weight, a.data.id, b.vertex, b.data.weight, b.data.id);
            });
        }
    }

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
     **/
    void setNumberOfSelfLoops(Unsafe, count loops) { storedNumberOfSelfLoops = loops; }

    /* NODE MODIFIERS */

    /**
     * Add a new node to the graph and return it.
     * @return The new node.
     */
    node addNode();

    /**
     * Add numberOfNewNodes new nodes.
     * @param  numberOfNewNodes Number of new nodes.
     * @return The index of the last node added.
     */
    node addNodes(count numberOfNewNodes);

    /** Removing a vertex is not supported in dhb. This will not remove the node but all neighbors
     * of the given node @a v
     * TODO:  In the future, we might want to make removeNode to use the concurrent version of
     * forNeighbors. But in that case, we need to make this function thread safe.
     * Incoming as well as outgoing edges will be removed.
     *
     * @param u Node.
     */
    void removeNode(node v);

    /**
     * Check if node @a v exists in the graph.
     *
     * @param v Node.
     * @return @c true if @a v exists, @c false otherwise.
     */

    bool hasNode(node v) const noexcept;

    /** NODE PROPERTIES **/
    /**
     * Returns the number of outgoing neighbors of @a v.
     *
     * @param v Node.
     * @return The number of outgoing neighbors.
     * @note The existence of the node is not checked. Calling this function with a non-existing
     * node results in a segmentation fault. Node existence can be checked by calling hasNode(u).
     */
    count degree(node v) const;

    /**
     * Get the number of incoming neighbors of @a v.
     *
     * @param v Node.
     * @return The number of incoming neighbors.
     * @note If the graph is not directed, the outgoing degree is returned.
     * @note The existence of the node is not checked. Calling this function with a non-existing
     * node results in a segmentation fault. Node existence can be checked by calling hasNode(u).
     */
    count degreeIn(node v) const;

    /**
     * Get the number of outgoing neighbors of @a v.
     *
     * @param v Node.
     * @return The number of outgoing neighbors.
     * @note The existence of the node is not checked. Calling this function with a non-existing
     * node results in a segmentation fault. Node existence can be checked by calling hasNode(u).
     */
    count degreeOut(node v) const { return degree(v); }

    /**
     * Check whether @a v is isolated, i.e. degree is 0.
     * @param v Node. If graph is directed, the function runs in O(n) since it needs to iterate over
     * all edges.
     * @return @c true if the node is isolated (= degree is 0)
     */
    bool isIsolated(node v) const;

    /**
     * Returns the weighted degree of @a u.
     *
     * @param u Node.
     * @param countSelfLoopsTwice If set to true, self-loops will be counted twice.
     *
     * @return Weighted degree of @a u.
     */
    edgeweight weightedDegree(node u, bool countSelfLoopsTwice = false) const;

    /**
     * Returns the weighted in-degree of @a u.
     *
     * @param u Node.
     * @param countSelfLoopsTwice If set to true, self-loops will be counted twice.
     *
     * @return Weighted in-degree of @a v.
     */
    edgeweight weightedDegreeIn(node u, bool countSelfLoopsTwice = false) const;

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
     * @param weight Optional edge weight.
     * @return @c true if edge has been added, false otherwise (in case checkMultiEdge is set to
     * true and the new edge would have been a multi-edge.)
     */
    bool addEdge(node u, node v, edgeweight ew = defaultEdgeWeight);

    /**
     * @brief Adds a collection of weighted edges to the DHBGraph.
     *
     * This function takes a vector of `WeightedEdge` objects and attempts to add them to the
     * `DHBGraph` instance. It processes the edges in parallel. Depending on the value of
     * `do_update`, the function either updates existing edges with new data or inserts new edges.
     * After adding the edges, edge ids of the new edges will not be updated. Please call
     * indexEdges(true) to index all the edges if needed.
     *
     * @param[in] weighted_edges A vector of `WeightedEdge` objects to be added to the graph. This
     * parameter is an rvalue reference.
     * @param[in] do_update A boolean flag indicating whether existing edges should be updated. If
     * `true`, the function will update existing edges; if `false`, it will only insert new edges.
     *
     * TODO: Document how to check for successful insertions.
     */
    void addEdges(std::vector<WeightedEdge> &&edges, bool do_update,
                  unsigned int num_threads = std::thread::hardware_concurrency());

    /**
     * @brief Converts a vector of `Edge` objects to `WeightedEdge` and adds them to the `DHBGraph`.
     *
     * @param[in] edges A vector of `Edge` objects to be converted and added to the graph. This
     * parameter is an rvalue reference.
     * @param[in] do_update A boolean flag indicating whether existing edges should be updated. If
     * `true`, the function will update existing edges; if `false`, it will only insert new edges.
     *
     * @return `true` if all edges were successfully added, `false` otherwise.
     */
    void addEdges(std::vector<Edge> &&edges, bool do_update,
                  unsigned int num_threads = std::thread::hardware_concurrency());

    /**
     *
     * Removes the undirected edge {@a u,@a v}.
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * DHB: Does not support maintaining compact edges or maintain sorted edges.
     */
    bool removeEdge(node u, node v);

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
    std::pair<count, count> removeAdjacentEdges(node u, Condition condition, bool edgesIn = false);

    /**
     * Removes all self-loops in the graph.
     */
    void removeSelfLoops();

    /**
     * Swaps edge target of edge (@a source_a, @a target_a) with edge target of
     * edge {@a source_b, @a target_b} such that (@a source_a, @a target_b) and
     * (@a source_b, @a target_a).
     *
     * Edge weights or edge ids are preserved.
     *
     * A std::runtime_error() is thrown if:
     * - edge (@a source_a, @a target_a) does not exist
     * - edge (@a source_b, @a target_b) does not exist
     * - edge (@a source_a, @a target_b) does already exist
     * - edge (@a source_b, @a target_a) does already exist
     */
    void swapEdge(node source_a, node target_a, node source_b, node target_b);

    /**
     * Checks if undirected edge {@a u,@a v} exists in the graph.
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @return <code>true</code> if the edge exists, <code>false</code>
     * otherwise.
     */
    bool hasEdge(node u, node v) const noexcept;

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
    bool isEmpty() const noexcept { return 0u == m_dhb_graph.vertices_count(); }

    /**
     * Return the number of nodes in the graph.
     * @return The number of nodes.
     */
    count numberOfNodes() const noexcept;

    /**
     * Return the number of edges in the graph.
     * @return The number of edges.
     */
    count numberOfEdges() const noexcept;

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
    index upperNodeIdBound() const noexcept;

    /**
     * Return edge weight of edge {@a u,@a v}. Returns 0 if edge does not
     * exist.
     * BEWARE: Running time is \Theta(1) for DHB implementation but \Theta(deg(u)) for fallback
     * NetworKit graph DS.
     *
     * @param u Endpoint of edge.
     * @param v Endpoint of edge.
     * @return Edge weight of edge {@a u,@a v} or 0 if edge does not exist.
     */
    edgeweight weight(node u, node v) const;

    /**
     * Set the weight of an edge. If the edge does not exist,
     * it will be inserted.
     *
     * @param[in]	u	endpoint of edge
     * @param[in]	v	endpoint of edge
     * @param[in]	weight	edge weight
     */
    void setWeight(node u, node v, edgeweight ew);

    /**
     * Set the weight to the i-th neighbour of u.
     *
     * @param[in]	u	endpoint of edge
     * @param[in]	i	index of the nexight
     * @param[in]	weight	edge weight
     */
    void setWeightAtIthNeighbor(Unsafe, node u, index i, edgeweight ew);

    /**
     * Increase the weight of an edge. If the edge does not exist,
     * it will be inserted.
     *
     * @param[in]	u	endpoint of edge
     * @param[in]	v	endpoint of edge
     * @param[in]	weight	edge weight
     */
    void increaseWeight(node u, node v, edgeweight ew);

    /* SUMS */

    /**
     * Returns the sum of all edge weights.
     * @return The sum of all edge weights.
     */
    edgeweight totalEdgeWeight() const noexcept;

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
     *
     * TODO: Write more test for edgeIterator and nodeIterator.
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
    NeighborRange<false> neighborRange(node u) const {
        assert(u < upperNodeIdBound());
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
    NeighborWeightRange<false> weightNeighborRange(node u) const {
        assert(isWeighted());
        assert(u < upperNodeIdBound());
        return NeighborWeightRange<false>(*this, u);
    }

    /**
     * Default Graph class: Returns the index of node v in the array of outgoing edges of node u.
     * DHB: Returns the index of node v in the neighbors of outgoing edges of node u.
     * @param u Node
     * @param v Node
     * @return index of node v in the array of outgoing edges of node u.
     * TODO: We might want to support concurrency for this function in the future.
     */
    index indexOfNeighbor(node u, node v) const;

    /**
     * Return the i-th (outgoing) neighbor of @a u.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return @a i-th (outgoing) neighbor of @a u, or @c none if no such
     * neighbor exists.
     */
    node getIthNeighbor(node u, index i) const;

    /**
     * Return the weight to the i-th (outgoing) neighbor of @a u.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return @a edge weight to the i-th (outgoing) neighbor of @a u, or @c +inf if no such
     * neighbor exists.
     */
    edgeweight getIthNeighborWeight(node u, index i) const;

    /**
     * Get i-th (outgoing) neighbor of @a u and the corresponding edge weight.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return pair: i-th (outgoing) neighbor of @a u and the corresponding
     * edge weight, or @c defaultEdgeWeight if unweighted.
     */
    std::pair<node, edgeweight> getIthNeighborWithWeight(node u, index i) const;

    /**
     * Get i-th (outgoing) neighbor of @a u and the corresponding edge weight.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return pair: i-th (outgoing) neighbor of @a u and the corresponding
     * edge weight, or @c defaultEdgeWeight if unweighted.
     *
     * This function is not supported for DHB. Please use getIthNeighborWithWeight(node u, index i).
     */
    std::pair<node, edgeweight> getIthNeighborWithWeight(Unsafe, node u, index i) const;

    /**
     * Get i-th (outgoing) neighbor of @a u and the corresponding edge id.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return pair: i-th (outgoing) neighbor of @a u and the corresponding
     * edge id, or @c none if no such neighbor exists.
     */
    std::pair<node, edgeid> getIthNeighborWithId(node u, index i) const;

    /* NODE ITERATORS */

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

    /** Iterate over all nodes of the graph and call @a handle (lambda
     * closure) as long as @a condition remains true. This allows for breaking
     * from a node loop.
     * TODO: Check if there's a way to use omp parallel for in the loop without introducing
     * atomic variable.
     * TODO: We may want to have explict parallel in the function name.
     * @param condition Returning <code>false</code> breaks the loop.
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename C, typename L>
    void forNodesWhile(C condition, L handle) const;

    /**
     * Iterate randomly over all nodes of the graph and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L>
    void forNodesInRandomOrder(L handle) const;

    /**
     * Iterate over all undirected pairs of nodes and call @a handle (lambda
     * closure).
     *
     * @param handle Takes parameters <code>(node, node)</code>.
     */
    template <typename L>
    void forNodePairs(L handle) const;

    /**
     * Iterate over all undirected pairs of nodes in parallel and call @a
     * handle (lambda closure).
     *
     * @param handle Takes parameters <code>(node, node)</code>.
     */
    template <typename L>
    void forNodePairsParallel(L handle) const;

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
    void forNeighborsOf(node u, L handle) const;

    /**
     * Iterate over all neighbors of a node and call @a handle (lamdba
     * closure) concurrently.
     *
     * @param u Node.
     * @param handle Takes parameter <code>(node)</code> or <code>(node,
     * edgeweight)</code> which is a neighbor of @a u.
     * @note For directed graphs only outgoing edges from @a u are considered.
     * A node is its own neighbor if there is a self-loop.
     *
     */
    template <typename L>
    void forNeighborsOfParallel(node u, L handle) const;

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
    void forEdgesOf(node u, L handle) const;

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
    void forEdgesOfParallel(node u, L handle) const;
    /**
     * Iterate over all neighbors of a node and call handler (lamdba closure).
     * For directed graphs only incoming edges from u are considered.
     */
    template <typename L>
    void forInNeighborsOf(node u, L handle) const;

    /**
     * Parallel version of `forInNeighborsOf`
     */
    template <typename L>
    void forInNeighborsOfParallel(node u, L handle) const;

    /**
     * Iterate over all incoming edges of a node and call handler (lamdba
     * closure).
     * @note For undirected graphs all edges incident to u are also incoming
     * edges.
     *
     * Handle takes parameters (u, v) or (u, v, w) where w is the edge weight.
     */
    template <typename L>
    void forInEdgesOf(node u, L handle) const;

    /**
     * Parallel version of `forInEdgesOf`
     */
    template <typename L>
    void forInEdgesOfParallel(node u, L handle) const;

    /* REDUCTION ITERATORS */

    /**
     * Iterate in parallel over all nodes and sum (reduce +) the values
     * returned by the handler
     */
    template <typename L>
    double sumForNodesParallel(L handle) const;

    /**
     * Iterate in parallel over all edges and sum (reduce +) the values
     * returned by the handler
     */
    template <typename L>
    double sumForEdgesParallel(L handle) const;
};

/* NODE ITERATORS */

template <typename L>
void DHBGraph::forNodes(L handle) const {
    for (node v = 0; v < m_dhb_graph.vertices_count(); ++v) {
        handle(v);
    }
}

template <typename L>
void DHBGraph::parallelForNodes(L handle) const {
#pragma omp parallel for schedule(guided)
    for (omp_index v = 0; v < static_cast<omp_index>(m_dhb_graph.vertices_count()); ++v) {
        handle(static_cast<node>(v));
    }
}

template <typename C, typename L>
void DHBGraph::forNodesWhile(C condition, L handle) const {
    for (dhb::Vertex v = 0; v < m_dhb_graph.vertices_count(); ++v) {
        if (!condition()) {
            break;
        }
        handle(v);
    }
}

template <typename L>
void DHBGraph::forNodesInRandomOrder(L handle) const {
    std::vector<node> randVec;
    randVec.reserve(numberOfNodes());
    forNodes([&](node u) { randVec.push_back(u); });
    std::shuffle(randVec.begin(), randVec.end(), Aux::Random::getURNG());
    for (node v : randVec) {
        handle(v);
    }
}

template <typename L>
void DHBGraph::balancedParallelForNodes(L handle) const {
// TODO: define min block size (and test it!)
#pragma omp parallel for schedule(guided)
    for (omp_index v = 0; v < static_cast<omp_index>(m_dhb_graph.vertices_count()); ++v) {
        handle(v);
    }
}

template <typename L>
void DHBGraph::forNodePairs(L handle) const {
    for (node u = 0; u < m_dhb_graph.vertices_count(); ++u) {
        for (node v = u + 1; v < m_dhb_graph.vertices_count(); ++v) {
            handle(u, v);
        }
    }
}

template <typename L>
void DHBGraph::forNodePairsParallel(L handle) const {
#pragma omp parallel for schedule(guided)
    for (omp_index u = 0; u < static_cast<omp_index>(m_dhb_graph.vertices_count()); ++u) {
        auto neighborhood = m_dhb_graph.neighbors(u);
        for (auto target = neighborhood.begin(); target != neighborhood.end(); ++target) {
            handle(u, target->vertex());
        }
    }
}

/* EDGE ITERATORS */

// /* HELPERS */

// implementation for graphIsDirected == true
template <bool graphIsDirected>
inline bool DHBGraph::useEdgeInIteration(node /* u */, node /* v */) const {
    return true;
}

template <bool hasWeights, bool graphHasEdgeIds>
std::tuple<edgeweight, edgeid> DHBGraph::getDHBEdgeData(node u, node v) const {
    auto const w = hasWeights ? weight(u, v) : defaultEdgeWeight;
    auto const id = graphHasEdgeIds ? edgeId(u, v) : none;
    return std::make_tuple(w, id);
}

// implementation for graphIsDirected == false
template <>
inline bool DHBGraph::useEdgeInIteration<false>(node u, node v) const {
    return u >= v;
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
void DHBGraph::forOutEdgesOfImpl(node u, L handle) const {
    assert(u < upperNodeIdBound());
    auto neighbors = m_dhb_graph.neighbors(u);

    for (size_t i = 0; i < neighbors.degree(); ++i) {
        dhb::Vertex const v = neighbors[i].vertex();
        auto const [w, id] = getDHBEdgeData<hasWeights, graphHasEdgeIds>(u, v);
        edgeLambda<L>(handle, u, v, w, id);
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
void DHBGraph::forOutEdgesOfImplParallel(node u, L handle) const {
    assert(u < upperNodeIdBound());
    auto neighbors = m_dhb_graph.neighbors(u);

#pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(neighbors.degree()); ++i) {
        dhb::Vertex const v = neighbors[i].vertex();
        auto const [w, id] = getDHBEdgeData<hasWeights, graphHasEdgeIds>(u, v);
        edgeLambda<L>(handle, u, v, w, id);
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
void DHBGraph::forInEdgesOfImpl(node u, L handle) const {
    assert(u < upperNodeIdBound());
    if (graphIsDirected) {
        for (dhb::Vertex from = 0; from < m_dhb_graph.vertices_count(); ++from) {
            for (auto it = neighborRange(from).begin(); it != neighborRange(from).end(); ++it) {
                dhb::Vertex const to = *it;
                if (u == to) {
                    auto const [w, id] = getDHBEdgeData<hasWeights, graphHasEdgeIds>(from, to);
                    edgeLambda<L>(handle, to, from, w, id);
                }
            }
        }
    } else {
        forNeighborsOf(u, handle);
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline void DHBGraph::forInEdgesOfImplParallel(node u, L handle) const {

    assert(u < upperNodeIdBound());
    if (graphIsDirected) {
#pragma omp parallel for schedule(guided)
        for (omp_index from = 0; from < static_cast<omp_index>(m_dhb_graph.vertices_count());
             ++from) {
            for (auto it = neighborRange(from).begin(); it != neighborRange(from).end(); ++it) {
                dhb::Vertex const to = *it;
                if (u == to) {
                    auto const [w, id] = getDHBEdgeData<hasWeights, graphHasEdgeIds>(from, to);
                    edgeLambda<L>(handle, to, from, w, id);
                }
            }
        }
    } else {
        forNeighborsOfParallel(u, handle);
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
void DHBGraph::forEdgeImpl(L handle) const {
    for (dhb::Vertex u = 0; u < m_dhb_graph.vertices_count(); ++u) {
        auto neighbors = m_dhb_graph.neighbors(u);
        for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
            dhb::Vertex const v = it->vertex();
            // If the graph is undirected, this function only call the handle on edges where u > v.
            // This avoids processing each edge twice (once for each direction).
            // However, in some cases, we might want to iterate over all edges
            // in the undirected graph, regardless of direction.
            if (useEdgeInIteration<graphIsDirected>(u, v)) {
                auto const [w, id] = getDHBEdgeData<hasWeights, graphHasEdgeIds>(u, v);
                edgeLambda<L>(handle, u, v, w, id);
            }
        }
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
void DHBGraph::forEdgesImplParallel(L handle) const {
#pragma omp parallel for schedule(guided)
    for (omp_index u = 0; u < static_cast<omp_index>(m_dhb_graph.vertices_count()); ++u) {
        auto neighbors = m_dhb_graph.neighbors(u);
        for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
            dhb::Vertex const v = it->vertex();
            // if the graph is not directed, only call the handle on edges that u > v.
            // So we don't want to use this sometimes, when we want to iterate all edges in an
            // undirected graph
            if (useEdgeInIteration<graphIsDirected>(u, v)) {
                auto const [w, id] = getDHBEdgeData<hasWeights, graphHasEdgeIds>(u, v);
                edgeLambda<L>(handle, u, v, w, id);
            }
        }
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
double DHBGraph::parallelSumForEdgesImpl(L handle) const {
    double total_sum = 0.0;
#pragma omp parallel for schedule(guided) reduction(+ : total_sum)
    for (omp_index u = 0; u < static_cast<omp_index>(m_dhb_graph.vertices_count()); ++u) {
        auto neighbors = m_dhb_graph.neighbors(u);

        for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
            dhb::Vertex const v = it->vertex();
            auto const [w, id] = getDHBEdgeData<hasWeights, graphHasEdgeIds>(u, v);
            if (u == v && !directed) {
                total_sum += 2 * edgeLambda<L>(handle, u, v, w, id);
            } else {
                total_sum += edgeLambda<L>(handle, u, v, w, id);
            }
        }
    }
    if (!directed) {
        total_sum /= 2;
    }
    return total_sum;
}

template <typename L>
void DHBGraph::forEdges(L handle) const {
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

template <typename L>
void DHBGraph::parallelForEdges(L handle) const {
    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edgeIds
        forEdgesImplParallel<false, false, false, L>(handle);
        break;

    case 1: // weighted,   undirected, no edgeIds
        forEdgesImplParallel<false, true, false, L>(handle);
        break;

    case 2: // unweighted, directed, no edgeIds
        forEdgesImplParallel<true, false, false, L>(handle);
        break;

    case 3: // weighted, directed, no edgeIds
        forEdgesImplParallel<true, true, false, L>(handle);
        break;

    case 4: // unweighted, undirected, with edgeIds
        forEdgesImplParallel<false, false, true, L>(handle);
        break;

    case 5: // weighted,   undirected, with edgeIds
        forEdgesImplParallel<false, true, true, L>(handle);
        break;

    case 6: // unweighted, directed, with edgeIds
        forEdgesImplParallel<true, false, true, L>(handle);
        break;

    case 7: // weighted,   directed, with edgeIds
        forEdgesImplParallel<true, true, true, L>(handle);
        break;
    }
}

/* NEIGHBORHOOD ITERATORS */

template <typename L>
void DHBGraph::forNeighborsOf(node u, L handle) const {
    forEdgesOf(u, handle);
}

template <typename L>
void DHBGraph::forEdgesOf(node u, L handle) const {
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

template <typename L>
void DHBGraph::forNeighborsOfParallel(node u, L handle) const {
    forEdgesOfParallel(u, handle);
}

template <typename L>
void DHBGraph::forEdgesOfParallel(node u, L handle) const {
    switch (weighted + 2 * edgesIndexed) {
    case 0: // not weighted, no edge ids
        forOutEdgesOfImplParallel<true, false, false, L>(u, handle);
        break;

    case 1: // weighted, no edge ids
        forOutEdgesOfImplParallel<true, true, false, L>(u, handle);
        break;

    case 2: // not weighted, with edge ids
        forOutEdgesOfImplParallel<true, false, true, L>(u, handle);
        break;

    case 3: // weighted, with edge ids
        forOutEdgesOfImplParallel<true, true, true, L>(u, handle);
        break;
    }
}

template <typename L>
void DHBGraph::forInNeighborsOf(node u, L handle) const {
    forInEdgesOf(u, handle);
}

template <typename L>
void DHBGraph::forInNeighborsOfParallel(node u, L handle) const {
    forInEdgesOfParallel(u, handle);
}

template <typename L>
void DHBGraph::forInEdgesOfParallel(node u, L handle) const {
    switch (weighted + 2 * directed + 4 * edgesIndexed) {
    case 0: // unweighted, undirected, no edge ids
        forInEdgesOfImplParallel<false, false, false, L>(u, handle);
        break;

    case 1: // weighted, undirected, no edge ids
        forInEdgesOfImplParallel<false, true, false, L>(u, handle);
        break;

    case 2: // unweighted, directed, no edge ids
        forInEdgesOfImplParallel<true, false, false, L>(u, handle);
        break;

    case 3: // weighted, directed, no edge ids
        forInEdgesOfImplParallel<true, true, false, L>(u, handle);
        break;

    case 4: // unweighted, undirected, with edge ids
        forInEdgesOfImplParallel<false, false, true, L>(u, handle);
        break;

    case 5: // weighted, undirected, with edge ids
        forInEdgesOfImplParallel<false, true, true, L>(u, handle);
        break;

    case 6: // unweighted, directed, with edge ids
        forInEdgesOfImplParallel<true, false, true, L>(u, handle);
        break;

    case 7: // weighted, directed, with edge ids
        forInEdgesOfImplParallel<true, true, true, L>(u, handle);
        break;
    }
}

template <typename L>
void DHBGraph::forInEdgesOf(node u, L handle) const {
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

template <typename L>
double DHBGraph::sumForNodesParallel(L handle) const {
    double sum = 0.0;

#pragma omp parallel for schedule(guided) reduction(+ : sum)
    for (omp_index v = 0; v < static_cast<omp_index>(m_dhb_graph.vertices_count()); ++v) {
        sum += handle(v);
    }

    return sum;
}

template <typename L>
double DHBGraph::sumForEdgesParallel(L handle) const {
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

template <typename Condition>
std::pair<count, count> DHBGraph::removeAdjacentEdges(node u, Condition condition, bool edgesIn) {
    std::vector<std::pair<dhb::Vertex, dhb::Vertex>> edgesToRemove;
    count removedEdges = 0;
    count removedSelfLoops = 0;

    if (!edgesIn) {
        forNeighborsOf(u, [&](node v) {
            if (condition(v)) {
                edgesToRemove.emplace_back(u, v);
            }
        });
    } else {
        forInNeighborsOf(u, [&](node v) {
            if (condition(v)) {
                edgesToRemove.emplace_back(v, u);
            }
        });
    }

    // actual removal
    for (auto const &edge : edgesToRemove) {
        bool const remove_success = removeEdge(edge.first, edge.second);
        bool const is_loop = edge.first == edge.second;
        if (remove_success) {
            removedSelfLoops += is_loop;
            removedEdges += !is_loop;
        }
    }
    if (!isDirected()) {
        return {removedEdges >> 1, removedSelfLoops};
    }
    return {removedEdges, removedSelfLoops};
}

} /* namespace NetworKit */

#endif // NETWORKIT_GRAPH_DHB_GRAPH_HPP_
