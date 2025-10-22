/*
 * EdgeIterators.hpp
 *
 *  Created on: 16.09.2024
 */

#ifndef NETWORKIT_GRAPH_EDGE_ITERATORS_HPP_
#define NETWORKIT_GRAPH_EDGE_ITERATORS_HPP_

#include <type_traits>

#include <networkit/Globals.hpp>
#include <networkit/graph/DynamicGraphUtils.hpp>
#include <networkit/graph/NodeIterators.hpp>

namespace NetworKit {

template <template <class, class> class GraphType, class NodeType, class EdgeWeightType>
class EdgeIteratorBase {
protected:
    const GraphType<NodeType, EdgeWeightType> *G;
    NodeIteratorBase<GraphType, NodeType, EdgeWeightType> nodeIter;
    index i{none};

public:
    EdgeIteratorBase(const GraphType<NodeType, EdgeWeightType> *G,
                     NodeIteratorBase<GraphType, NodeType, EdgeWeightType> nodeIter)
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

    // A valid edge might be defined differently for different graph types
    bool validEdge() const noexcept {
        return G->isDirected() || (*nodeIter <= G->getIthNeighbor(Unsafe{}, *nodeIter, i));
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

    bool operator==(const EdgeIteratorBase &rhs) const noexcept {
        return nodeIter == rhs.nodeIter && i == rhs.i;
    }

    bool operator!=(const EdgeIteratorBase &rhs) const noexcept { return !(*this == rhs); }
};

template <class NodeType>
struct is_weighted_edge_specialization : std::false_type {};

template <class NodeType, class EdgeWeightType>
struct is_weighted_edge_specialization<WeightedEdge<NodeType, EdgeWeightType>> : std::true_type {};

/**
 * Class to iterate over the edges of the given graph type. If the graph is undirected, operator*()
 * returns the edges (u, v) s.t. u <= v.
 */
template <template <class, class> class GraphType, class NodeType, class EdgeWeightType,
          class IterEdgeWeightType>
class EdgeWeightTypeIterator : public EdgeIteratorBase<GraphType, NodeType, EdgeWeightType> {
    using EdgeIteratorBase<GraphType, NodeType, EdgeWeightType>::nodeIter;
    using EdgeIteratorBase<GraphType, NodeType, EdgeWeightType>::G;
    using EdgeIteratorBase<GraphType, NodeType, EdgeWeightType>::i;

    Edge<NodeType> starOperator(std::false_type) const {
        return Edge<NodeType>(*nodeIter, G->getIthNeighbor(Unsafe{}, *nodeIter, i));
    }

    WeightedEdge<NodeType, EdgeWeightType> starOperator(std::true_type) const {
        return WeightedEdge<NodeType, EdgeWeightType>(
            *nodeIter, G->getIthNeighbor(Unsafe{}, *nodeIter, i),
            G->getIthNeighborWeight(Unsafe{}, *nodeIter, i));
    }

public:
    // The value type of the edges (i.e. a pair). Returned by operator*().
    using value_type = IterEdgeWeightType;

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
    using self = EdgeWeightTypeIterator;

    EdgeWeightTypeIterator(const GraphType<NodeType, EdgeWeightType> *G,
                           NodeIteratorBase<GraphType, NodeType, EdgeWeightType> nodeIter)
        : EdgeIteratorBase<GraphType, NodeType, EdgeWeightType>(G, nodeIter) {}

    EdgeWeightTypeIterator() : EdgeIteratorBase<GraphType, NodeType, EdgeWeightType>() {}

    bool operator==(const EdgeWeightTypeIterator &rhs) const noexcept {
        return this->EdgeIteratorBase<GraphType, NodeType, EdgeWeightType>::operator==(
            static_cast<EdgeIteratorBase<GraphType, NodeType, EdgeWeightType>>(rhs));
    }

    bool operator!=(const EdgeWeightTypeIterator &rhs) const noexcept { return !(*this == rhs); }

    // The returned edge depends on both the type of edge and graph
    IterEdgeWeightType operator*() const noexcept {
        return this->starOperator(is_weighted_edge_specialization<IterEdgeWeightType>{});
    }

    EdgeWeightTypeIterator &operator++() {
        EdgeIteratorBase<GraphType, NodeType, EdgeWeightType>::nextEdge();
        return *this;
    }

    EdgeWeightTypeIterator operator++(int) {
        const auto tmp = *this;
        ++(*this);
        return tmp;
    }

    EdgeWeightTypeIterator operator--() {
        EdgeIteratorBase<GraphType, NodeType, EdgeWeightType>::prevEdge();
        return *this;
    }

    EdgeWeightTypeIterator operator--(int) {
        const auto tmp = *this;
        --(*this);
        return tmp;
    }
};

/**
 * Wrapper class to iterate over a range of the edges.
 */
template <template <class, class> class GraphType, class NodeType, class EdgeWeightType,
          class IterEdgeWeightType>
class EdgeWeightTypeRange {

    const GraphType<NodeType, EdgeWeightType> *G;

public:
    EdgeWeightTypeRange(const GraphType<NodeType, EdgeWeightType> &G) : G(&G) {}

    EdgeWeightTypeRange() : G(nullptr){};

    ~EdgeWeightTypeRange() = default;

    EdgeWeightTypeIterator<GraphType, NodeType, EdgeWeightType, IterEdgeWeightType> begin() const {
        assert(G);
        return EdgeWeightTypeIterator<GraphType, NodeType, EdgeWeightType, IterEdgeWeightType>(
            G, G->nodeRange().begin());
    }

    EdgeWeightTypeIterator<GraphType, NodeType, EdgeWeightType, IterEdgeWeightType> end() const {
        assert(G);
        return EdgeWeightTypeIterator<GraphType, NodeType, EdgeWeightType, IterEdgeWeightType>(
            G, G->nodeRange().end());
    }
};

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_EDGE_ITERATORS_HPP_
