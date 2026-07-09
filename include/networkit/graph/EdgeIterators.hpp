/*
 * EdgeIterators.hpp
 *
 *  Created on: 16.09.2024
 */

#ifndef NETWORKIT_GRAPH_EDGE_ITERATORS_HPP_
#define NETWORKIT_GRAPH_EDGE_ITERATORS_HPP_

#include <type_traits>

#include <networkit/Globals.hpp>
#include <networkit/graph/EdgeUtils.hpp>
#include <networkit/graph/NodeIterators.hpp>

namespace NetworKit {

template <class GraphType, class NodeT, class EdgeWeightT>
class EdgeIteratorBase {
protected:
    const GraphType *G;
    NodeIteratorBase<GraphType, NodeT, EdgeWeightT> nodeIter;
    index i{none};

public:
    EdgeIteratorBase(const GraphType *G, NodeIteratorBase<GraphType, NodeT, EdgeWeightT> nodeIter)
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

template <class T>
struct is_weighted_edge_specialization : std::false_type {};

template <class NodeT, class EdgeWeightT>
struct is_weighted_edge_specialization<WeightedEdgeT<NodeT, EdgeWeightT>> : std::true_type {};

/**
 * Class to iterate over the edges of the given graph type. If the graph is undirected, operator*()
 * returns the edges (u, v) s.t. u <= v.
 */
template <class GraphType, class NodeT, class EdgeWeightT, class IterEdgeWeightT>
class EdgeWeightTIterator : public EdgeIteratorBase<GraphType, NodeT, EdgeWeightT> {
    using EdgeIteratorBase<GraphType, NodeT, EdgeWeightT>::nodeIter;
    using EdgeIteratorBase<GraphType, NodeT, EdgeWeightT>::G;
    using EdgeIteratorBase<GraphType, NodeT, EdgeWeightT>::i;

    EdgeT<NodeT> starOperator(std::false_type) const {
        return EdgeT<NodeT>(*nodeIter, G->getIthNeighbor(Unsafe{}, *nodeIter, i));
    }

    WeightedEdgeT<NodeT, EdgeWeightT> starOperator(std::true_type) const {
        return WeightedEdgeT<NodeT, EdgeWeightT>(*nodeIter,
                                                 G->getIthNeighbor(Unsafe{}, *nodeIter, i),
                                                 G->getIthNeighborWeight(Unsafe{}, *nodeIter, i));
    }

public:
    // The value type of the edges (i.e. a pair). Returned by operator*().
    using value_type = IterEdgeWeightT;

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
    using self = EdgeWeightTIterator;

    EdgeWeightTIterator(const GraphType *G,
                        NodeIteratorBase<GraphType, NodeT, EdgeWeightT> nodeIter)
        : EdgeIteratorBase<GraphType, NodeT, EdgeWeightT>(G, nodeIter) {}

    EdgeWeightTIterator() : EdgeIteratorBase<GraphType, NodeT, EdgeWeightT>() {}

    bool operator==(const EdgeWeightTIterator &rhs) const noexcept {
        return this->EdgeIteratorBase<GraphType, NodeT, EdgeWeightT>::operator==(
            static_cast<EdgeIteratorBase<GraphType, NodeT, EdgeWeightT>>(rhs));
    }

    bool operator!=(const EdgeWeightTIterator &rhs) const noexcept { return !(*this == rhs); }

    // The returned edge depends on both the type of edge and graph
    IterEdgeWeightT operator*() const noexcept {
        return this->starOperator(is_weighted_edge_specialization<IterEdgeWeightT>{});
    }

    EdgeWeightTIterator &operator++() {
        EdgeIteratorBase<GraphType, NodeT, EdgeWeightT>::nextEdge();
        return *this;
    }

    EdgeWeightTIterator operator++(int) {
        const auto tmp = *this;
        ++(*this);
        return tmp;
    }

    EdgeWeightTIterator operator--() {
        EdgeIteratorBase<GraphType, NodeT, EdgeWeightT>::prevEdge();
        return *this;
    }

    EdgeWeightTIterator operator--(int) {
        const auto tmp = *this;
        --(*this);
        return tmp;
    }
};

/**
 * Wrapper class to iterate over a range of the edges.
 */
template <class GraphType, class NodeT, class EdgeWeightT, class IterEdgeWeightT>
class EdgeWeightTRange {

    const GraphType *G;

public:
    EdgeWeightTRange(const GraphType &G) : G(&G) {}

    EdgeWeightTRange() : G(nullptr) {};

    ~EdgeWeightTRange() = default;

    EdgeWeightTIterator<GraphType, NodeT, EdgeWeightT, IterEdgeWeightT> begin() const {
        assert(G);
        return EdgeWeightTIterator<GraphType, NodeT, EdgeWeightT, IterEdgeWeightT>(
            G, G->nodeRange().begin());
    }

    EdgeWeightTIterator<GraphType, NodeT, EdgeWeightT, IterEdgeWeightT> end() const {
        assert(G);
        return EdgeWeightTIterator<GraphType, NodeT, EdgeWeightT, IterEdgeWeightT>(
            G, G->nodeRange().end());
    }
};

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_EDGE_ITERATORS_HPP_
