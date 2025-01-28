/*
 * NodeIterators.hpp
 *
 *  Created on: 16.09.2024
 */

#ifndef NETWORKIT_GRAPH_NODE_ITERATORS_HPP_
#define NETWORKIT_GRAPH_NODE_ITERATORS_HPP_

#include <cassert>
#include <iterator>

#include <networkit/Globals.hpp>

namespace NetworKit {

/**
 * Class to iterate over the nodes of a graph.
 */
template <typename GraphType>
class NodeIteratorBase {

    const GraphType *G;
    node u{none};

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
    using self = NodeIteratorBase;

    NodeIteratorBase(const GraphType *G, node u) : G(G), u(u) {
        if (!G->hasNode(u) && u < G->upperNodeIdBound()) {
            ++(*this);
        }
    }

    /**
     * @brief WARNING: This constructor is required for Python and should not be used as the
     * iterator is not initialized.
     */
    NodeIteratorBase() : G(nullptr) {}

    ~NodeIteratorBase() = default;

    NodeIteratorBase &operator++() {
        assert(u < G->upperNodeIdBound());
        do {
            ++u;
        } while (!(G->hasNode(u) || u >= G->upperNodeIdBound()));
        return *this;
    }

    NodeIteratorBase operator++(int) {
        const auto tmp = *this;
        ++(*this);
        return tmp;
    }

    NodeIteratorBase operator--() {
        assert(u);
        do {
            --u;
        } while (!G->hasNode(u));
        return *this;
    }

    NodeIteratorBase operator--(int) {
        const auto tmp = *this;
        --(*this);
        return tmp;
    }

    bool operator==(const NodeIteratorBase &rhs) const noexcept { return u == rhs.u; }

    bool operator!=(const NodeIteratorBase &rhs) const noexcept { return !(*this == rhs); }

    node operator*() const noexcept {
        assert(u < G->upperNodeIdBound());
        return u;
    }
};

/**
 * Wrapper class to iterate over a range of the nodes of a graph.
 */
template <typename GraphType>
class NodeRangeBase {

    const GraphType *G;

public:
    NodeRangeBase(const GraphType &G) : G(&G) {}

    NodeRangeBase() : G(nullptr){};

    ~NodeRangeBase() = default;

    NodeIteratorBase<GraphType> begin() const noexcept {
        assert(G);
        return NodeIteratorBase(G, node{0});
    }

    NodeIteratorBase<GraphType> end() const noexcept {
        assert(G);
        return NodeIteratorBase(G, G->upperNodeIdBound());
    }
};

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_NODE_ITERATORS_HPP_
