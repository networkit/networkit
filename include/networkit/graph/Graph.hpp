/*
 * Graph.hpp
 *
 *  Created on: 01.06.2014
 *      Author: Christian Staudt
 *              Klara Reichard <klara.reichard@gmail.com>
 *              Marvin Ritter <marvin.ritter@gmail.com>
 */

#ifndef NETWORKIT_GRAPH_GRAPH_HPP_
#define NETWORKIT_GRAPH_GRAPH_HPP_

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <omp.h>
#include <queue>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <typeindex>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/ArrayTools.hpp>
#include <networkit/auxiliary/FunctionTraits.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph_helper/Edge.hpp>

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
class Graph final {

    // graph attributes
    //!< current number of nodes
    count n;
    //!< current number of edges
    count m;

    //!< current number of self loops, edges which have the same origin and
    //!< target
    count storedNumberOfSelfLoops;

    //!< current upper bound of node ids, z will be the id of the next node
    node z;
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
    std::vector<std::vector<node>> inEdges;
    //!< (outgoing) edges, for each edge (u, v) v is saved in outEdges[u] and
    //!< for undirected also u in outEdges[v]
    std::vector<std::vector<node>> outEdges;

    //!< only used for directed graphs, same schema as inEdges
    std::vector<std::vector<edgeweight>> inEdgeWeights;
    //!< same schema (and same order!) as outEdges
    std::vector<std::vector<edgeweight>> outEdgeWeights;

    //!< only used for directed graphs, same schema as inEdges
    std::vector<std::vector<edgeid>> inEdgeIds;
    //!< same schema (and same order!) as outEdges
    std::vector<std::vector<edgeid>> outEdgeIds;

private:
    // base class for all node (and edge) attribute
    // storages with attribute type info
    // independent of the attribute type, holds bookkeeping info only:
    // - attribute name
    // - type info of derived (real storage holding) classes
    // - which indices are valid
    // - number of valid indices
    // - the associated graph (who knows, which nodes/edges exist)
    // - the validity of the whole storage (initially true, false after detach)
    // all indexed accesses by NetworKit::index: synonym both for node and edgeid

    class PerNode {
    public:
        static constexpr bool edges = false;
    };
    class PerEdge {
    public:
        static constexpr bool edges = true;
    };

    template <typename NodeOrEdge>
    class AttributeStorageBase { // alias ASB
    public:
        AttributeStorageBase(const Graph *graph, std::string name, std::type_index type)
            : name{std::move(name)}, type{type}, theGraph{graph}, validStorage{true} {
            checkPremise(); // node for PerNode, theGraph.hasEdgeIds() for PerEdges
        }

        void invalidateStorage() { validStorage = false; }

        const std::string &getName() const noexcept { return name; }

        std::type_index getType() const noexcept { return type; }

        bool isValid(index n) const noexcept { return n < valid.size() && valid[n]; }

        // Called by Graph when node/edgeid n is deleted.
        void invalidate(index n) {
            if (isValid(n)) {
                valid[n] = false;
                --validElements;
            }
        }

    protected:
        void markValid(index n) {
            indexOK(n); // specialized for node/edgeid
            if (n >= valid.size())
                valid.resize(n + 1);
            if (!valid[n]) {
                valid[n] = true;
                ++validElements;
            }
        }

        void checkIndex(index n) const {
            indexOK(n);
            if (!isValid(n)) {
                throw std::runtime_error("Invalid attribute value");
            }
        }

    private:
        std::string name;
        std::type_index type;
        std::vector<bool> valid; // For each node/edgeid: whether attribute is set or not.

    protected:
        void indexOK(index n) const;
        void checkPremise() const;
        index validElements = 0;
        const Graph *theGraph;
        bool validStorage; // Validity of the whole storage

    }; // class AttributeStorageBase

    template <typename NodeOrEdge>
    using ASB = AttributeStorageBase<NodeOrEdge>;

    template <typename NodeOrEdge, typename T, bool isConst>
    class Attribute;

    template <typename NodeOrEdge, template <typename> class Base, typename T>
    class AttributeStorage : public Base<NodeOrEdge> {
    public:
        AttributeStorage(const Graph *theGraph, std::string name)
            : Base<NodeOrEdge>{theGraph, std::move(name), typeid(T)} {}

        void resize(index i) {
            if (i >= values.size())
                values.resize(i + 1);
        }

        auto size() const noexcept { return this->validElements; }

        void set(index i, T &&v) {
            this->markValid(i);
            resize(i);
            values[i] = std::move(v);
        }

        // instead of returning an std::optional (C++17) we provide these
        // C++14 options
        // (1) throw an exception when invalid:
        T get(index i) const { // may throw
            this->checkIndex(i);
            return values[i];
        }

        // (2) give default value when invalid:
        T get(index i, T defaultT) const noexcept {
            if (i >= values.size() || !this->isValid(i))
                return defaultT;
            return values[i];
        }

        friend Attribute<NodeOrEdge, T, true>;
        friend Attribute<NodeOrEdge, T, false>;

    private:
        using Base<NodeOrEdge>::theGraph;
        std::vector<T> values; // the real attribute storage
    };                         // class AttributeStorage<NodeOrEdge, Base, T>

    template <typename NodeOrEdge, typename T, bool isConst>
    class Attribute {
    public:
        using AttributeStorage_type =
            std::conditional_t<isConst, const AttributeStorage<NodeOrEdge, ASB, T>,
                               AttributeStorage<NodeOrEdge, ASB, T>>;
        class Iterator {
        public:
            // The value type of the attribute. Returned by
            // operator*().
            using value_type = T;

            // Reference to the value_type, required by STL.
            using reference = std::conditional_t<isConst, const value_type &, value_type &>;

            // Pointer to the value_type, required by STL.
            using pointer = std::conditional_t<isConst, const value_type *, value_type *>;

            // STL iterator category.
            using iterator_category = std::forward_iterator_tag;

            // Signed integer type of the result of subtracting two pointers,
            // required by STL.
            using difference_type = ptrdiff_t;

            Iterator() : storage{nullptr}, idx{0} {}
            Iterator(AttributeStorage_type *storage) : storage{storage}, idx{0} {
                if (storage) {
                    nextValid();
                }
            }

            Iterator &nextValid() {
                while (storage && !storage->isValid(idx)) {
                    if (idx >= storage->values.size()) {
                        storage = nullptr;
                        return *this;
                    }
                    ++idx;
                }
                return *this;
            }

            Iterator &operator++() {
                if (!storage) {
                    throw std::runtime_error("Invalid attribute iterator");
                }
                ++idx;
                return nextValid();
            }

            auto operator*() const {
                if (!storage) {
                    throw std::runtime_error("Invalid attribute iterator");
                }
                return std::make_pair(idx, storage->values[idx]);
            }

            bool operator==(Iterator const &iter) const noexcept {
                if (storage == nullptr && iter.storage == nullptr) {
                    return true;
                }
                return storage == iter.storage && idx == iter.idx;
            }

            bool operator!=(Iterator const &iter) const noexcept { return !(*this == iter); }

        private:
            AttributeStorage_type *storage;
            index idx;
        }; // class Iterator

    private:
        class IndexProxy {
            // a helper class for distinguished read and write on an indexed
            // attribute
            // operator[] on an attribute yields an IndexProxy holding
            // location and index of access
            //    - casting an IndexProxy to the attribute type reads the value
            //    - assigning to it (operator=) writes the value
        public:
            IndexProxy(AttributeStorage_type *storage, index idx) : storage{storage}, idx{idx} {}

            // reading at idx
            operator T() const {
                storage->checkIndex(idx);
                return storage->values[idx];
            }

            // writing at idx
            template <bool ic = isConst>
            std::enable_if_t<!ic, T> &operator=(T &&other) {
                storage->set(idx, std::move(other));
                return storage->values[idx];
            }

        private:
            AttributeStorage_type *storage;
            index idx;
        }; // class IndexProxy
    public:
        explicit Attribute(std::shared_ptr<AttributeStorage_type> ownedStorage = nullptr)
            : ownedStorage{ownedStorage}, valid{ownedStorage != nullptr} {}

        Attribute(Attribute const &other) : ownedStorage{other.ownedStorage}, valid{other.valid} {}

        template <bool ic = isConst, std::enable_if_t<ic, int> = 0>
        Attribute(Attribute<NodeOrEdge, T, false> const &other)
            : ownedStorage{other.ownedStorage}, valid{other.valid} {}

        Attribute &operator=(Attribute other) {
            this->swap(other);
            return *this;
        }

        void swap(Attribute &other) {
            std::swap(ownedStorage, other.ownedStorage);
            std::swap(valid, other.valid);
        }

        Attribute(Attribute &&other) noexcept
            : ownedStorage{std::move(other.ownedStorage)}, valid{other.valid} {
            other.valid = false;
        }

        template <bool ic = isConst, std::enable_if_t<ic, int> = 0>
        Attribute(Attribute<NodeOrEdge, T, false> &&other) noexcept
            : ownedStorage{std::move(other.ownedStorage)}, valid{other.valid} {
            other.valid = false;
        }

        auto begin() const {
            checkAttribute();
            return Iterator(ownedStorage.get()).nextValid();
        }

        auto end() const { return Iterator(nullptr); }

        auto size() const noexcept { return ownedStorage->size(); }

        template <bool ic = isConst>
        std::enable_if_t<!ic> set(index i, T v) {
            checkAttribute();
            ownedStorage->set(i, std::move(v));
        }

        template <bool ic = isConst>
        std::enable_if_t<!ic> set2(node u, node v, T t) {
            static_assert(NodeOrEdge::edges, "attribute(u,v) for edges only");
            set(ownedStorage->theGraph->edgeId(u, v), t);
        }

        auto get(index i) const {
            checkAttribute();
            return ownedStorage->get(i);
        }

        auto get2(node u, node v) const {
            static_assert(NodeOrEdge::edges, "attribute(u,v) for edges only");
            return get(ownedStorage->theGraph->edgeId(u, v));
        }

        auto get(index i, T defaultT) const {
            checkAttribute();
            return ownedStorage->get(i, defaultT);
        }

        auto get2(node u, node v, T defaultT) const {
            static_assert(NodeOrEdge::edges, "attribute(u,v) for edges only");
            return get(ownedStorage->theGraph->edgeId(u, v), defaultT);
        }

        IndexProxy operator[](index i) const {
            checkAttribute();
            return IndexProxy(ownedStorage.get(), i);
        }

        IndexProxy operator()(node u, node v) const {
            static_assert(NodeOrEdge::edges, "attribute(u,v) for edges only");
            checkAttribute();
            return IndexProxy(ownedStorage.get(), ownedStorage->theGraph->edgeId(u, v));
        }

        void checkAttribute() const {
            if (!ownedStorage->validStorage)
                throw std::runtime_error("Invalid attribute");
        }

        auto getName() const {
            checkAttribute();
            return ownedStorage->getName();
        }

        void write(std::string const &filename) const {
            std::ofstream out(filename);
            if (!out)
                ERROR("cannot open ", filename, " for writing");

            for (auto it = begin(); it != end(); ++it) {
                auto pair = *it;
                auto n = pair.first;  // node/edgeid
                auto v = pair.second; // value
                out << n << "\t" << v << "\n";
            }
            out.close();
        }

        template <bool ic = isConst>
        std::enable_if_t<!ic> read(const std::string &filename) {
            std::ifstream in(filename);
            if (!in) {
                ERROR("cannot open ", filename, " for reading");
            }
            index n; // node/edgeid
            T v;     // value
            std::string line;
            while (std::getline(in, line)) {
                std::istringstream istring(line);
                if constexpr (std::is_same_v<T, std::string>) {
                    istring >> n >> std::ws;
                    std::getline(istring, v);
                } else {
                    istring >> n >> v;
                }
                set(n, v);
            }
        }

    private:
        std::shared_ptr<AttributeStorage_type> ownedStorage;
        bool valid;
    }; // class Attribute

    template <typename NodeOrEdge>
    class AttributeMap {
        friend Graph;
        const Graph *theGraph;

    public:
        std::unordered_map<std::string, std::shared_ptr<ASB<NodeOrEdge>>> attrMap;

        AttributeMap(const Graph *g) : theGraph{g} {}

        auto find(std::string const &name) {
            auto it = attrMap.find(name);
            if (it == attrMap.end()) {
                throw std::runtime_error("No such attribute");
            }
            return it;
        }

        auto find(std::string const &name) const {
            auto it = attrMap.find(name);
            if (it == attrMap.end()) {
                throw std::runtime_error("No such attribute");
            }
            return it;
        }

        template <typename T>
        auto attach(const std::string &name) {
            auto ownedPtr =
                std::make_shared<AttributeStorage<NodeOrEdge, ASB, T>>(theGraph, std::string{name});
            auto insertResult = attrMap.emplace(ownedPtr->getName(), ownedPtr);
            auto success = insertResult.second;
            if (!success) {
                throw std::runtime_error("Attribute with same name already exists");
            }
            return Attribute<NodeOrEdge, T, false>{ownedPtr};
        }

        void detach(const std::string &name) {
            auto it = find(name);
            auto storage = it->second.get();
            storage->invalidateStorage();
            it->second.reset();
            attrMap.erase(name);
        }

        template <typename T>
        auto get(const std::string &name) {
            auto it = find(name);
            if (it->second.get()->getType() != typeid(T))
                throw std::runtime_error("Type mismatch in Attributes().get()");
            return Attribute<NodeOrEdge, T, false>{
                std::static_pointer_cast<AttributeStorage<NodeOrEdge, ASB, T>>(it->second)};
        }

        template <typename T>
        auto get(const std::string &name) const {
            auto it = find(name);
            if (it->second.get()->getType() != typeid(T))
                throw std::runtime_error("Type mismatch in Attributes().get()");
            return Attribute<NodeOrEdge, T, true>{
                std::static_pointer_cast<const AttributeStorage<NodeOrEdge, ASB, T>>(it->second)};
        }

    }; // class AttributeMap

    AttributeMap<PerNode> nodeAttributeMap;
    AttributeMap<PerEdge> edgeAttributeMap;

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

    using NodeIntAttribute = Attribute<PerNode, int, false>;
    using NodeDoubleAttribute = Attribute<PerNode, double, false>;
    using NodeStringAttribute = Attribute<PerNode, std::string, false>;

    using EdgeIntAttribute = Attribute<PerEdge, int, false>;
    using EdgeDoubleAttribute = Attribute<PerEdge, double, false>;
    using EdgeStringAttribute = Attribute<PerEdge, std::string, false>;

private:
    /**
     * Returns the index of node u in the array of incoming edges of node v.
     * (for directed graphs inEdges is searched, while for indirected outEdges
     * is searched, which gives the same result as indexInOutEdgeArray).
     */
    index indexInInEdgeArray(node v, node u) const;

    /**
     * Returns the index of node v in the array of outgoing edges of node u.
     */
    index indexInOutEdgeArray(node u, node v) const;

    /**
     * Computes the weighted in/out degree of node @a u.
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
     * Returns the edge weight of the outgoing edge of index i in the outgoing
     * edges of node u
     * @param u The node
     * @param i The index
     * @return The weight of the outgoing edge or defaultEdgeWeight if the graph
     * is unweighted
     */
    template <bool hasWeights>
    inline edgeweight getOutEdgeWeight(node u, index i) const;

    /**
     * Returns the edge weight of the incoming edge of index i in the incoming
     * edges of node u
     *
     * @param u The node
     * @param i The index in the incoming edge array
     * @return The weight of the incoming edge
     */
    template <bool hasWeights>
    inline edgeweight getInEdgeWeight(node u, index i) const;

    /**
     * Returns the edge id of the edge of index i in the outgoing edges of node
     * u
     *
     * @param u The node
     * @param i The index in the outgoing edges
     * @return The edge id
     */
    template <bool graphHasEdgeIds>
    inline edgeid getOutEdgeId(node u, index i) const;

    /**
     * Returns the edge id of the edge of index i in the incoming edges of node
     * u
     *
     * @param u The node
     * @param i The index in the incoming edges of u
     * @return The edge id
     */
    template <bool graphHasEdgeIds>
    inline edgeid getInEdgeId(node u, index i) const;

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
    inline void forOutEdgesOfImpl(node u, L handle) const;

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
    inline void forInEdgesOfImpl(node u, L handle) const;

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
    /**
     * Class to iterate over the nodes of a graph.
     */
    class NodeIterator {

        const Graph *G;
        node u;

    public:
        // The value type of the nodes (i.e. nodes). Returned by
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
        using self = NodeIterator;

        NodeIterator(const Graph *G, node u) : G(G), u(u) {
            if (!G->hasNode(u) && u < G->upperNodeIdBound()) {
                ++(*this);
            }
        }

        /**
         * @brief WARNING: This constructor is required for Python and should not be used as the
         * iterator is not initialized.
         */
        NodeIterator() : G(nullptr) {}

        ~NodeIterator() = default;

        NodeIterator &operator++() {
            assert(u < G->upperNodeIdBound());
            do {
                ++u;
            } while (!(G->hasNode(u) || u >= G->upperNodeIdBound()));
            return *this;
        }

        NodeIterator operator++(int) {
            const auto tmp = *this;
            ++(*this);
            return tmp;
        }

        NodeIterator operator--() {
            assert(u);
            do {
                --u;
            } while (!G->hasNode(u));
            return *this;
        }

        NodeIterator operator--(int) {
            const auto tmp = *this;
            --(*this);
            return tmp;
        }

        bool operator==(const NodeIterator &rhs) const noexcept { return u == rhs.u; }

        bool operator!=(const NodeIterator &rhs) const noexcept { return !(*this == rhs); }

        node operator*() const noexcept {
            assert(u < G->upperNodeIdBound());
            return u;
        }
    };

    /**
     * Wrapper class to iterate over a range of the nodes of a graph.
     */
    class NodeRange {

        const Graph *G;

    public:
        NodeRange(const Graph &G) : G(&G) {}

        NodeRange() : G(nullptr){};

        ~NodeRange() = default;

        NodeIterator begin() const noexcept {
            assert(G);
            return NodeIterator(G, node{0});
        }

        NodeIterator end() const noexcept {
            assert(G);
            return NodeIterator(G, G->upperNodeIdBound());
        }
    };

    // Necessary for friendship with EdgeIteratorBase.
    class EdgeIterator;
    class EdgeWeightIterator;

    class EdgeIteratorBase {
        friend class EdgeIterator;
        friend class EdgeWeightIterator;

        const Graph *G;
        NodeIterator nodeIter;
        index i;

        bool validEdge() const noexcept {
            return G->isDirected() || (*nodeIter <= G->outEdges[*nodeIter][i]);
        }

        void nextEdge() {
            do {
                if (++i >= G->degree(*nodeIter)) {
                    i = 0;
                    do {
                        assert(nodeIter != G->nodeRange().end());
                        ++nodeIter;
                        if (nodeIter == G->nodeRange().end()) {
                            return;
                        }
                    } while (!G->degree(*nodeIter));
                }
            } while (!validEdge());
        }

        void prevEdge() {
            do {
                if (!i) {
                    do {
                        assert(nodeIter != G->nodeRange().begin());
                        --nodeIter;
                    } while (!G->degree(*nodeIter));

                    i = G->degree(*nodeIter);
                }
                --i;
            } while (!validEdge());
        }

        EdgeIteratorBase(const Graph *G, NodeIterator nodeIter)
            : G(G), nodeIter(nodeIter), i(index{0}) {
            if (nodeIter != G->nodeRange().end() && !G->degree(*nodeIter)) {
                nextEdge();
            }
        }

        /**
         * @brief WARNING: This constructor is required for Python and should not be used as the
         * iterator is not initialized.
         */
        EdgeIteratorBase() : G(nullptr) {}

        virtual ~EdgeIteratorBase() = default;

        bool operator==(const EdgeIteratorBase &rhs) const noexcept {
            return nodeIter == rhs.nodeIter && i == rhs.i;
        }

        bool operator!=(const EdgeIteratorBase &rhs) const noexcept { return !(*this == rhs); }
    };

    /**
     * Class to iterate over the edges of the graph. If the graph is undirected, operator*()
     * returns the edges (u, v) s.t. u <= v.
     */
    class EdgeIterator : public EdgeIteratorBase {

    public:
        // The value type of the edges (i.e. a pair). Returned by operator*().
        using value_type = Edge;

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
        using self = EdgeIterator;

        EdgeIterator(const Graph *G, NodeIterator nodeIter) : EdgeIteratorBase(G, nodeIter) {}

        EdgeIterator() : EdgeIteratorBase() {}

        bool operator==(const EdgeIterator &rhs) const noexcept {
            return this->EdgeIteratorBase::operator==(static_cast<EdgeIteratorBase>(rhs));
        }

        bool operator!=(const EdgeIterator &rhs) const noexcept { return !(*this == rhs); }

        Edge operator*() const noexcept {
            assert(nodeIter != G->nodeRange().end());
            return Edge(*nodeIter, G->outEdges[*nodeIter][i]);
        }

        EdgeIterator &operator++() {
            nextEdge();
            return *this;
        }

        EdgeIterator operator++(int) {
            const auto tmp = *this;
            ++(*this);
            return tmp;
        }

        EdgeIterator operator--() {
            prevEdge();
            return *this;
        }

        EdgeIterator operator--(int) {
            const auto tmp = *this;
            --(*this);
            return tmp;
        }
    };

    /**
     * Class to iterate over the edges of the graph and their weights. If the graph is undirected,
     * operator*() returns a WeightedEdge struct with u <= v.
     */
    class EdgeWeightIterator : public EdgeIteratorBase {
    public:
        // The value type of the edges and their weights (i.e. WeightedEdge). Returned by
        // operator*().
        using value_type = WeightedEdge;

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
        using self = EdgeWeightIterator;

        EdgeWeightIterator(const Graph *G, NodeIterator nodeIter) : EdgeIteratorBase(G, nodeIter) {}

        /**
         * @brief WARNING: This constructor is required for Python and should not be used as the
         * iterator is not initialized.
         */
        EdgeWeightIterator() : EdgeIteratorBase() {}

        bool operator==(const EdgeWeightIterator &rhs) const noexcept {
            return this->EdgeIteratorBase::operator==(static_cast<EdgeIteratorBase>(rhs));
        }

        bool operator!=(const EdgeWeightIterator &rhs) const noexcept { return !(*this == rhs); }

        EdgeWeightIterator &operator++() {
            nextEdge();
            return *this;
        }

        EdgeWeightIterator operator++(int) {
            const auto tmp = *this;
            ++(*this);
            return tmp;
        }

        EdgeWeightIterator operator--() {
            prevEdge();
            return *this;
        }

        EdgeWeightIterator operator--(int) {
            const auto tmp = *this;
            --(*this);
            return tmp;
        }

        WeightedEdge operator*() const noexcept {
            assert(nodeIter != G->nodeRange().end());
            return WeightedEdge(*nodeIter, G->outEdges[*nodeIter][i],
                                G->isWeighted() ? G->outEdgeWeights[*nodeIter][i] : 1);
        }
    };

    /**
     * Wrapper class to iterate over a range of the edges of a graph.
     */
    class EdgeRange {

        const Graph *G;

    public:
        EdgeRange(const Graph &G) : G(&G) {}

        EdgeRange() : G(nullptr){};

        ~EdgeRange() = default;

        EdgeIterator begin() const {
            assert(G);
            return EdgeIterator(G, G->nodeRange().begin());
        }

        EdgeIterator end() const {
            assert(G);
            return EdgeIterator(G, G->nodeRange().end());
        }
    };

    /**
     * Wrapper class to iterate over a range of the edges of a graph and their weights.
     */
    class EdgeWeightRange {

        const Graph *G;

    public:
        EdgeWeightRange(const Graph &G) : G(&G) {}

        EdgeWeightRange() : G(nullptr){};

        ~EdgeWeightRange() = default;

        EdgeWeightIterator begin() const {
            assert(G);
            return EdgeWeightIterator(G, G->nodeRange().begin());
        }

        EdgeWeightIterator end() const {
            assert(G);
            return EdgeWeightIterator(G, G->nodeRange().end());
        }
    };

    /**
     * Class to iterate over the in/out neighbors of a node.
     */
    class NeighborIterator {

        std::vector<node>::const_iterator nIter;

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

        NeighborIterator(std::vector<node>::const_iterator nodesIter) : nIter(nodesIter) {}

        /**
         * @brief WARNING: This contructor is required for Python and should not be used as the
         * iterator is not initialized.
         */
        NeighborIterator() {}

        NeighborIterator &operator++() {
            ++nIter;
            return *this;
        }

        NeighborIterator operator++(int) {
            const auto tmp = *this;
            ++nIter;
            return tmp;
        }

        NeighborIterator operator--() {
            const auto tmp = *this;
            --nIter;
            return tmp;
        }

        NeighborIterator operator--(int) {
            --nIter;
            return *this;
        }

        bool operator==(const NeighborIterator &rhs) const { return nIter == rhs.nIter; }

        bool operator!=(const NeighborIterator &rhs) const { return !(nIter == rhs.nIter); }

        node operator*() const { return *nIter; }
    };

    /**
     * Class to iterate over the in/out neighbors of a node including the edge
     * weights. Values are std::pair<node, edgeweight>.
     */
    class NeighborWeightIterator {

        std::vector<node>::const_iterator nIter;
        std::vector<edgeweight>::const_iterator wIter;

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

        NeighborWeightIterator(std::vector<node>::const_iterator nodesIter,
                               std::vector<edgeweight>::const_iterator weightIter)
            : nIter(nodesIter), wIter(weightIter) {}

        /**
         * @brief WARNING: This contructor is required for Python and should not be used as the
         * iterator is not initialized.
         */
        NeighborWeightIterator() {}

        NeighborWeightIterator &operator++() {
            ++nIter;
            ++wIter;
            return *this;
        }

        NeighborWeightIterator operator++(int) {
            const auto tmp = *this;
            ++(*this);
            return tmp;
        }

        NeighborWeightIterator operator--() {
            --nIter;
            --wIter;
            return *this;
        }

        NeighborWeightIterator operator--(int) {
            const auto tmp = *this;
            --(*this);
            return tmp;
        }

        bool operator==(const NeighborWeightIterator &rhs) const {
            return nIter == rhs.nIter && wIter == rhs.wIter;
        }

        bool operator!=(const NeighborWeightIterator &rhs) const { return !(*this == rhs); }

        std::pair<node, edgeweight> operator*() const { return std::make_pair(*nIter, *wIter); }
    };

    /**
     * Wrapper class to iterate over a range of the neighbors of a node within
     * a for loop.
     */
    template <bool InEdges = false>
    class NeighborRange {
        const Graph *G;
        node u;

    public:
        NeighborRange(const Graph &G, node u) : G(&G), u(u) { assert(G.hasNode(u)); };

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
     * Values are std::pair<node, edgeweight>.
     */
    template <bool InEdges = false>
    class NeighborWeightRange {

        const Graph *G;
        node u;

    public:
        NeighborWeightRange(const Graph &G, node u) : G(&G), u(u) { assert(G.hasNode(u)); };

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
    Graph(count n = 0, bool weighted = false, bool directed = false, bool edgesIndexed = false);

    template <class EdgeMerger = std::plus<edgeweight>>
    Graph(const Graph &G, bool weighted, bool directed, bool edgesIndexed = false,
          EdgeMerger edgeMerger = std::plus<edgeweight>())
        : n(G.n), m(G.m), storedNumberOfSelfLoops(G.storedNumberOfSelfLoops), z(G.z),
          omega(edgesIndexed ? G.omega : 0), t(G.t), weighted(weighted), directed(directed),
          edgesIndexed(edgesIndexed), // edges are not indexed by default
          exists(G.exists),

          // let the following be empty for the start, we fill them later
          inEdges(0), outEdges(0), inEdgeWeights(0), outEdgeWeights(0), inEdgeIds(0), outEdgeIds(0),

          // empty node attribute map as last member for this graph
          nodeAttributeMap(this), edgeAttributeMap(this) {

        if (G.isDirected() == directed) {
            // G.inEdges might be empty (if G is undirected), but
            // that's fine
            inEdges = G.inEdges;
            outEdges = G.outEdges;

            // copy weights if needed
            if (weighted) {
                if (G.isWeighted()) {
                    // just copy from G, again either both graphs are directed or both are
                    // undirected
                    inEdgeWeights = G.inEdgeWeights;
                    outEdgeWeights = G.outEdgeWeights;
                } else {
                    // G has no weights, set defaultEdgeWeight for all edges
                    if (directed) {
                        inEdgeWeights.resize(z);
                        for (node u = 0; u < z; u++) {
                            inEdgeWeights[u].resize(G.inEdges[u].size(), defaultEdgeWeight);
                        }
                    }

                    outEdgeWeights.resize(z);
                    for (node u = 0; u < z; u++) {
                        outEdgeWeights[u].resize(outEdges[u].size(), defaultEdgeWeight);
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
            outEdges.resize(z);
            if (weighted)
                outEdgeWeights.resize(z);
            if (G.hasEdgeIds() && edgesIndexed)
                outEdgeIds.resize(z);
            G.balancedParallelForNodes([&](node u) {
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
                        outEdgeWeights[u].resize(outEdges[u].size(), defaultEdgeWeight);
                    }
                }
                if (G.hasEdgeIds() && edgesIndexed) {
                    // copy both out and in edges ids into our new outEdgesIds
                    outEdgeIds[u].reserve(G.outEdgeIds[u].size() + G.inEdgeIds[u].size());
                    outEdgeIds[u].insert(outEdgeIds[u].end(), G.outEdgeIds[u].begin(),
                                         G.outEdgeIds[u].end());
                }
            });
            G.balancedParallelForNodes([&](node u) {
                // this is necessary to avoid multi edges, because both u -> v and v -> u can exist
                // in G
                count edgeSurplus = 0;
                for (count i = 0; i < G.inEdges[u].size(); ++i) {
                    node v = G.inEdges[u][i];
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
                                        : edgeMerger(defaultEdgeWeight, defaultEdgeWeight);
                            if (G.hasEdgeIds() && edgesIndexed)
                                outEdgeIds[u][j] = std::min(G.inEdgeIds[u][i], G.outEdgeIds[u][j]);
                        }
                        break;
                    }
                    if (!alreadyPresent) { // an equivalent out edge wasn't present so we add it
                        outEdges[u].push_back(v);
                        if (weighted)
                            outEdgeWeights[u].push_back(G.isWeighted() ? G.inEdgeWeights[u][i]
                                                                       : defaultEdgeWeight);
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
            inEdges = G.outEdges;
            outEdges = G.outEdges;
            if (weighted) {
                if (G.isWeighted()) {
                    inEdgeWeights = G.outEdgeWeights;
                    outEdgeWeights = G.outEdgeWeights;
                } else {
                    // initialize both inEdgeWeights and outEdgeWeights with the
                    // defaultEdgeWeight
                    inEdgeWeights.resize(z);
                    for (node u = 0; u < z; ++u) {
                        inEdgeWeights[u].resize(inEdges[u].size(), defaultEdgeWeight);
                    }
                    outEdgeWeights.resize(z);
                    for (node u = 0; u < z; ++u) {
                        outEdgeWeights[u].resize(outEdges[u].size(), defaultEdgeWeight);
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
    Graph(std::initializer_list<WeightedEdge> edges);

    /**
     * Create a graph as copy of @a other.
     * @param other The graph to copy.
     */
    Graph(const Graph &other) = default;

    /** Default move constructor */
    Graph(Graph &&other) noexcept = default;

    /** Default destructor */
    ~Graph() = default;

    /** Default move assignment operator */
    Graph &operator=(Graph &&other) noexcept = default;

    /** Default copy assignment operator */
    Graph &operator=(const Graph &other) = default;

    /**
     * Reserves memory in the node's edge containers for undirected graphs.
     *
     * @param u the node memory should be reserved for
     * @param size the amount of memory to reserve
     *
     * This function is thread-safe if called from different
     * threads on different nodes.
     */
    void preallocateUndirected(node u, size_t size);

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
    void preallocateDirected(node u, size_t outSize, size_t inSize);

    /**
     * Reserves memory in the node's edge containers for directed graphs.
     *
     * @param u the node memory should be reserved for
     * @param outSize the amount of memory to reserve for out edges
     *
     * This function is thread-safe if called from different
     * threads on different nodes.
     */
    void preallocateDirectedOutEdges(node u, size_t outSize);

    /**
     * Reserves memory in the node's edge containers for directed graphs.
     *
     * @param u the node memory should be reserved for
     * @param inSize the amount of memory to reserve for in edges
     *
     * This function is thread-safe if called from different
     * threads on different nodes.
     */
    void preallocateDirectedInEdges(node u, size_t inSize);

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
    node addNode();

    /**
     * Add numberOfNewNodes new nodes.
     * @param  numberOfNewNodes Number of new nodes.
     * @return The index of the last node added.
     */
    node addNodes(count numberOfNewNodes);

    /**
     * Remove a node @a v and all incident edges from the graph.
     *
     * Incoming as well as outgoing edges will be removed.
     *
     * @param u Node.
     */
    void removeNode(node v);

    /**
     * Removes out-going edges from node @u. If the graph is weighted and/or has edge ids, weights
     * and/or edge ids will also be removed.
     *
     * @param node u Node.
     */
    void removePartialOutEdges(Unsafe, node u) {
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
     * @param node u Node.
     */
    void removePartialInEdges(Unsafe, node u) {
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

    bool hasNode(node v) const noexcept { return (v < z) && this->exists[v]; }

    /**
     * Restores a previously deleted node @a v with its previous id in the
     * graph.
     *
     * @param v Node.
     *
     */

    void restoreNode(node v);

    /** NODE PROPERTIES **/
    /**
     * Returns the number of outgoing neighbors of @a v.
     *
     * @param v Node.
     * @return The number of outgoing neighbors.
     * @note The existence of the node is not checked. Calling this function with a non-existing
     * node results in a segmentation fault. Node existence can be checked by calling hasNode(u).
     */
    count degree(node v) const {
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
    count degreeIn(node v) const {
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
    count degreeOut(node v) const { return degree(v); }

    /**
     * Check whether @a v is isolated, i.e. degree is 0.
     * @param v Node.
     * @return @c true if the node is isolated (= degree is 0)
     */
    bool isIsolated(node v) const {
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
     * @param checkMultiEdge If true, this enables a check for a possible multi-edge.
     * @return @c true if edge has been added, false otherwise (in case checkMultiEdge is set to
     * true and the new edge would have been a multi-edge.)
     */
    bool addEdge(node u, node v, edgeweight ew = defaultEdgeWeight, bool checkMultiEdge = false);

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
     * @param weight Optional edge weight.
     * @param ew Optional edge weight.
     * @param index Optional edge index.
     * @param checkMultiEdge If true, this enables a check for a possible multi-edge.
     * @return @c true if edge has been added, false otherwise (in case checkMultiEdge is set to
     * true and the new edge would have been a multi-edge.)
     */
    bool addPartialEdge(Unsafe, node u, node v, edgeweight ew = defaultEdgeWeight,
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
     * @param checkMultiEdge If true, this enables a check for a possible multi-edge.
     * @return @c true if edge has been added, false otherwise (in case checkMultiEdge is set to
     * true and the new edge would have been a multi-edge.)
     */
    bool addPartialInEdge(Unsafe, node u, node v, edgeweight ew = defaultEdgeWeight,
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
     * @param checkMultiEdge If true, this enables a check for a possible multi-edge.
     * @return @c true if edge has been added, false otherwise (in case checkMultiEdge is set to
     * true and the new edge would have been a multi-edge.)
     */
    bool addPartialOutEdge(Unsafe, node u, node v, edgeweight ew = defaultEdgeWeight,
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
    void removeEdge(node u, node v);

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
    void swapEdge(node s1, node t1, node s2, node t2);

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
        t++;
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
     * Set the weight to the i-th incoming neighbour of u.
     *
     * @param[in]	u	endpoint of edge
     * @param[in]	i	index of the nexight
     * @param[in]	weight	edge weight
     */
    void setWeightAtIthInNeighbor(Unsafe, node u, index i, edgeweight ew);

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
    NeighborWeightRange<false> weightNeighborRange(node u) const {
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
    NeighborRange<true> inNeighborRange(node u) const {
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
    NeighborWeightRange<true> weightInNeighborRange(node u) const {
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
    index indexOfNeighbor(node u, node v) const { return indexInOutEdgeArray(u, v); }

    /**
     * Return the i-th (outgoing) neighbor of @a u.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return @a i-th (outgoing) neighbor of @a u, or @c none if no such
     * neighbor exists.
     */
    node getIthNeighbor(node u, index i) const {
        if (!hasNode(u) || i >= outEdges[u].size())
            return none;
        return outEdges[u][i];
    }

    /**
     * Return the i-th (incoming) neighbor of @a u.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeIn(u))
     * @return @a i-th (incoming) neighbor of @a u, or @c none if no such
     * neighbor exists.
     */
    node getIthInNeighbor(node u, index i) const {
        if (!hasNode(u) || i >= inEdges[u].size())
            return none;
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
    edgeweight getIthNeighborWeight(node u, index i) const {
        if (!hasNode(u) || i >= outEdges[u].size())
            return nullWeight;
        return isWeighted() ? outEdgeWeights[u][i] : defaultEdgeWeight;
    }

    /**
     * Get i-th (outgoing) neighbor of @a u and the corresponding edge weight.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return pair: i-th (outgoing) neighbor of @a u and the corresponding
     * edge weight, or @c defaultEdgeWeight if unweighted.
     */
    std::pair<node, edgeweight> getIthNeighborWithWeight(node u, index i) const {
        if (!hasNode(u) || i >= outEdges[u].size())
            return {none, none};
        return getIthNeighborWithWeight(unsafe, u, i);
    }

    /**
     * Get i-th (outgoing) neighbor of @a u and the corresponding edge weight.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return pair: i-th (outgoing) neighbor of @a u and the corresponding
     * edge weight, or @c defaultEdgeWeight if unweighted.
     */
    std::pair<node, edgeweight> getIthNeighborWithWeight(Unsafe, node u, index i) const {
        if (!isWeighted())
            return {outEdges[u][i], defaultEdgeWeight};
        return {outEdges[u][i], outEdgeWeights[u][i]};
    }

    /**
     * Get i-th (outgoing) neighbor of @a u and the corresponding edge id.
     *
     * @param u Node.
     * @param i index; should be in [0, degreeOut(u))
     * @return pair: i-th (outgoing) neighbor of @a u and the corresponding
     * edge id, or @c none if no such neighbor exists.
     */
    std::pair<node, edgeid> getIthNeighborWithId(node u, index i) const {
        assert(hasEdgeIds());
        if (!hasNode(u) || i >= outEdges[u].size())
            return {none, none};
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
     *
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
     * Iterate in parallel over all nodes of the graph and call handler
     * (lambda closure). Using schedule(guided) to remedy load-imbalances due
     * to e.g. unequal degree distribution.
     *
     * @param handle Takes parameter <code>(node)</code>.
     */
    template <typename L>
    void balancedParallelForNodes(L handle) const;

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
    void parallelForNodePairs(L handle) const;

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
     * Iterate over all neighbors of a node and call handler (lamdba closure).
     * For directed graphs only incoming edges from u are considered.
     */
    template <typename L>
    void forInNeighborsOf(node u, L handle) const;

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

    /* REDUCTION ITERATORS */

    /**
     * Iterate in parallel over all nodes and sum (reduce +) the values
     * returned by the handler
     */
    template <typename L>
    double parallelSumForNodes(L handle) const;

    /**
     * Iterate in parallel over all edges and sum (reduce +) the values
     * returned by the handler
     */
    template <typename L>
    double parallelSumForEdges(L handle) const;
};

/* NODE ITERATORS */

template <typename L>
void Graph::forNodes(L handle) const {
    for (node v = 0; v < z; ++v) {
        if (exists[v]) {
            handle(v);
        }
    }
}

template <typename L>
void Graph::parallelForNodes(L handle) const {
#pragma omp parallel for
    for (omp_index v = 0; v < static_cast<omp_index>(z); ++v) {
        if (exists[v]) {
            handle(v);
        }
    }
}

template <typename C, typename L>
void Graph::forNodesWhile(C condition, L handle) const {
    for (node v = 0; v < z; ++v) {
        if (exists[v]) {
            if (!condition()) {
                break;
            }
            handle(v);
        }
    }
}

template <typename L>
void Graph::forNodesInRandomOrder(L handle) const {
    std::vector<node> randVec;
    randVec.reserve(numberOfNodes());
    forNodes([&](node u) { randVec.push_back(u); });
    std::shuffle(randVec.begin(), randVec.end(), Aux::Random::getURNG());
    for (node v : randVec) {
        handle(v);
    }
}

template <typename L>
void Graph::balancedParallelForNodes(L handle) const {
// TODO: define min block size (and test it!)
#pragma omp parallel for schedule(guided)
    for (omp_index v = 0; v < static_cast<omp_index>(z); ++v) {
        if (exists[v]) {
            handle(v);
        }
    }
}

template <typename L>
void Graph::forNodePairs(L handle) const {
    for (node u = 0; u < z; ++u) {
        if (exists[u]) {
            for (node v = u + 1; v < z; ++v) {
                if (exists[v]) {
                    handle(u, v);
                }
            }
        }
    }
}

template <typename L>
void Graph::parallelForNodePairs(L handle) const {
#pragma omp parallel for schedule(guided)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        if (exists[u]) {
            for (node v = u + 1; v < z; ++v) {
                if (exists[v]) {
                    handle(u, v);
                }
            }
        }
    }
}

/* EDGE ITERATORS */

/* HELPERS */

template <typename T>
void erase(node u, index idx, std::vector<std::vector<T>> &vec);
// implementation for weighted == true
template <bool hasWeights>
inline edgeweight Graph::getOutEdgeWeight(node u, index i) const {
    return outEdgeWeights[u][i];
}

// implementation for weighted == false
template <>
inline edgeweight Graph::getOutEdgeWeight<false>(node, index) const {
    return defaultEdgeWeight;
}

// implementation for weighted == true
template <bool hasWeights>
inline edgeweight Graph::getInEdgeWeight(node u, index i) const {
    return inEdgeWeights[u][i];
}

// implementation for weighted == false
template <>
inline edgeweight Graph::getInEdgeWeight<false>(node, index) const {
    return defaultEdgeWeight;
}

// implementation for hasEdgeIds == true
template <bool graphHasEdgeIds>
inline edgeid Graph::getOutEdgeId(node u, index i) const {
    return outEdgeIds[u][i];
}

// implementation for hasEdgeIds == false
template <>
inline edgeid Graph::getOutEdgeId<false>(node, index) const {
    return none;
}

// implementation for hasEdgeIds == true
template <bool graphHasEdgeIds>
inline edgeid Graph::getInEdgeId(node u, index i) const {
    return inEdgeIds[u][i];
}

// implementation for hasEdgeIds == false
template <>
inline edgeid Graph::getInEdgeId<false>(node, index) const {
    return none;
}

// implementation for graphIsDirected == true
template <bool graphIsDirected>
inline bool Graph::useEdgeInIteration(node /* u */, node /* v */) const {
    return true;
}

// implementation for graphIsDirected == false
template <>
inline bool Graph::useEdgeInIteration<false>(node u, node v) const {
    return u >= v;
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline void Graph::forOutEdgesOfImpl(node u, L handle) const {
    for (index i = 0; i < outEdges[u].size(); ++i) {
        node v = outEdges[u][i];

        if (useEdgeInIteration<graphIsDirected>(u, v)) {
            edgeLambda<L>(handle, u, v, getOutEdgeWeight<hasWeights>(u, i),
                          getOutEdgeId<graphHasEdgeIds>(u, i));
        }
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline void Graph::forInEdgesOfImpl(node u, L handle) const {
    if (graphIsDirected) {
        for (index i = 0; i < inEdges[u].size(); i++) {
            node v = inEdges[u][i];

            if (useEdgeInIteration<true>(u, v)) {
                edgeLambda<L>(handle, u, v, getInEdgeWeight<hasWeights>(u, i),
                              getInEdgeId<graphHasEdgeIds>(u, i));
            }
        }
    } else {
        for (index i = 0; i < outEdges[u].size(); ++i) {
            node v = outEdges[u][i];

            if (useEdgeInIteration<true>(u, v)) {
                edgeLambda<L>(handle, u, v, getOutEdgeWeight<hasWeights>(u, i),
                              getOutEdgeId<graphHasEdgeIds>(u, i));
            }
        }
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline void Graph::forEdgeImpl(L handle) const {
    for (node u = 0; u < z; ++u) {
        forOutEdgesOfImpl<graphIsDirected, hasWeights, graphHasEdgeIds, L>(u, handle);
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline void Graph::parallelForEdgesImpl(L handle) const {
#pragma omp parallel for schedule(guided)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        forOutEdgesOfImpl<graphIsDirected, hasWeights, graphHasEdgeIds, L>(u, handle);
    }
}

template <bool graphIsDirected, bool hasWeights, bool graphHasEdgeIds, typename L>
inline double Graph::parallelSumForEdgesImpl(L handle) const {
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (omp_index u = 0; u < static_cast<omp_index>(z); ++u) {
        for (index i = 0; i < outEdges[u].size(); ++i) {
            node v = outEdges[u][i];

            // undirected, do not iterate over edges twice
            // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
            if (useEdgeInIteration<graphIsDirected>(u, v)) {
                sum += edgeLambda<L>(handle, u, v, getOutEdgeWeight<hasWeights>(u, i),
                                     getOutEdgeId<graphHasEdgeIds>(u, i));
            }
        }
    }

    return sum;
}

template <typename L>
void Graph::forEdges(L handle) const {
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
void Graph::parallelForEdges(L handle) const {
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

template <typename L>
void Graph::forNeighborsOf(node u, L handle) const {
    forEdgesOf(u, handle);
}

template <typename L>
void Graph::forEdgesOf(node u, L handle) const {
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
void Graph::forInNeighborsOf(node u, L handle) const {
    forInEdgesOf(u, handle);
}

template <typename L>
void Graph::forInEdgesOf(node u, L handle) const {
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
double Graph::parallelSumForNodes(L handle) const {
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (omp_index v = 0; v < static_cast<omp_index>(z); ++v) {
        if (exists[v]) {
            sum += handle(v);
        }
    }

    return sum;
}

template <typename L>
double Graph::parallelSumForEdges(L handle) const {
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
std::pair<count, count> Graph::removeAdjacentEdges(node u, Condition condition, bool edgesIn) {
    count removedEdges = 0;
    count removedSelfLoops = 0;

    // For directed graphs, this function is supposed to be called twice: one to remove out-edges,
    // and one to remove in-edges.
    auto &edges_ = edgesIn ? inEdges[u] : outEdges[u];
    for (index vi = 0; vi < edges_.size();) {
        if (condition(edges_[vi])) {
            const auto isSelfLoop = (edges_[vi] == u);
            removedSelfLoops += isSelfLoop;
            removedEdges += !isSelfLoop;
            edges_[vi] = edges_.back();
            edges_.pop_back();
            if (isWeighted()) {
                auto &weights_ = edgesIn ? inEdgeWeights[u] : outEdgeWeights[u];
                weights_[vi] = weights_.back();
                weights_.pop_back();
            }
            if (hasEdgeIds()) {
                auto &edgeIds_ = edgesIn ? inEdgeIds[u] : outEdgeIds[u];
                edgeIds_[vi] = edgeIds_.back();
                edgeIds_.pop_back();
            }
        } else {
            ++vi;
        }
    }

    return {removedEdges, removedSelfLoops};
}

template <class Lambda>
void Graph::sortEdges(Lambda lambda) {

    std::vector<std::vector<index>> indicesGlobal(omp_get_max_threads());

    const auto sortAdjacencyArrays = [&](node u, std::vector<node> &adjList,
                                         std::vector<edgeweight> &weights,
                                         std::vector<edgeid> &edgeIds) -> void {
        auto &indices = indicesGlobal[omp_get_thread_num()];
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
                return lambda(WeightedEdgeWithId{u, adjList[a], defaultEdgeWeight, edgeIds[a]},
                              WeightedEdgeWithId{u, adjList[b], defaultEdgeWeight, edgeIds[b]});
            });
        else
            std::sort(indices.begin(), indicesEnd, [&](auto a, auto b) -> bool {
                return lambda(WeightedEdgeWithId{u, adjList[a], defaultEdgeWeight, 0},
                              WeightedEdgeWithId{u, adjList[b], defaultEdgeWeight, 0});
            });

        Aux::ArrayTools::applyPermutation(adjList.begin(), adjList.end(), indices.begin());

        if (isWeighted())
            Aux::ArrayTools::applyPermutation(weights.begin(), weights.end(), indices.begin());

        if (hasEdgeIds())
            Aux::ArrayTools::applyPermutation(edgeIds.begin(), edgeIds.end(), indices.begin());
    };

    balancedParallelForNodes([&](const node u) {
        if (degree(u) < 2)
            return;

        std::vector<edgeweight> dummyEdgeWeights;
        std::vector<edgeid> dummyEdgeIds;
        sortAdjacencyArrays(u, outEdges[u], isWeighted() ? outEdgeWeights[u] : dummyEdgeWeights,
                            hasEdgeIds() ? outEdgeIds[u] : dummyEdgeIds);

        if (isDirected())
            sortAdjacencyArrays(u, inEdges[u], isWeighted() ? inEdgeWeights[u] : dummyEdgeWeights,
                                hasEdgeIds() ? inEdgeIds[u] : dummyEdgeIds);
    });
}

} /* namespace NetworKit */

#endif // NETWORKIT_GRAPH_GRAPH_HPP_
