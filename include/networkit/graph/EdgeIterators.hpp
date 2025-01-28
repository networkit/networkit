/*
 * EdgeIterators.hpp
 *
 *  Created on: 16.09.2024
 */

#ifndef NETWORKIT_GRAPH_EDGE_ITERATORS_HPP_
#define NETWORKIT_GRAPH_EDGE_ITERATORS_HPP_

#include <networkit/Globals.hpp>
#include <networkit/graph/NodeIterators.hpp>

namespace NetworKit {

template <typename GraphType>
class EdgeIteratorBase {

protected:
    const GraphType *G;
    NodeIteratorBase<GraphType> nodeIter;
    index i{none};

public:
    EdgeIteratorBase(const GraphType *G, NodeIteratorBase<GraphType> nodeIter)
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
    bool validEdge() const noexcept;

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

/**
 * Class to iterate over the edges of the given graph type. If the graph is undirected, operator*()
 * returns the edges (u, v) s.t. u <= v.
 */
template <typename GraphType, typename EdgeType>
class EdgeTypeIterator : public EdgeIteratorBase<GraphType> {

public:
    // The value type of the edges (i.e. a pair). Returned by operator*().
    using value_type = EdgeType;

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
    using self = EdgeTypeIterator;

    EdgeTypeIterator(const GraphType *G, NodeIteratorBase<GraphType> nodeIter)
        : EdgeIteratorBase<GraphType>(G, nodeIter) {}

    EdgeTypeIterator() : EdgeIteratorBase<GraphType>() {}

    bool operator==(const EdgeTypeIterator &rhs) const noexcept {
        return this->EdgeIteratorBase<GraphType>::operator==(
            static_cast<EdgeIteratorBase<GraphType>>(rhs));
    }

    bool operator!=(const EdgeTypeIterator &rhs) const noexcept { return !(*this == rhs); }

    // The returned edge depends on both the type of edge and graph
    EdgeType operator*() const noexcept;

    EdgeTypeIterator &operator++() {
        EdgeIteratorBase<GraphType>::nextEdge();
        return *this;
    }

    EdgeTypeIterator operator++(int) {
        const auto tmp = *this;
        ++(*this);
        return tmp;
    }

    EdgeTypeIterator operator--() {
        EdgeIteratorBase<GraphType>::prevEdge();
        return *this;
    }

    EdgeTypeIterator operator--(int) {
        const auto tmp = *this;
        --(*this);
        return tmp;
    }
};

/**
 * Wrapper class to iterate over a range of the edges.
 */
template <typename GraphType, typename EdgeType>
class EdgeTypeRange {

    const GraphType *G;

public:
    EdgeTypeRange(const GraphType &G) : G(&G) {}

    EdgeTypeRange() : G(nullptr){};

    ~EdgeTypeRange() = default;

    EdgeTypeIterator<GraphType, EdgeType> begin() const {
        assert(G);
        return EdgeTypeIterator<GraphType, EdgeType>(G, G->nodeRange().begin());
    }

    EdgeTypeIterator<GraphType, EdgeType> end() const {
        assert(G);
        return EdgeTypeIterator<GraphType, EdgeType>(G, G->nodeRange().end());
    }
};

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_EDGE_ITERATORS_HPP_
