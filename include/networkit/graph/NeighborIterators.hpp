/*
 * EdgeIterators.hpp
 *
 *  Created on: 14.05.2024
 */

#ifndef NETWORKIT_GRAPH_NEIGHBOR_ITERATORS_HPP_
#define NETWORKIT_GRAPH_NEIGHBOR_ITERATORS_HPP_

#include <networkit/Globals.hpp>
#include <networkit/graph/NodeIterators.hpp>

namespace NetworKit {

/**
 * Class to iterate over the in/out neighbors of a node.
 */
class NeighborIteratorBase {

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
    using self = NeighborIteratorBase;

    NeighborIteratorBase(std::vector<node>::const_iterator nodesIter) : nIter(nodesIter) {}

    /**
     * @brief WARNING: This contructor is required for Python and should not be used as the
     * iterator is not initialized.
     */
    NeighborIteratorBase() {}

    NeighborIteratorBase &operator++() {
        ++nIter;
        return *this;
    }

    NeighborIteratorBase operator++(int) {
        const auto tmp = *this;
        ++nIter;
        return tmp;
    }

    NeighborIteratorBase operator--() {
        const auto tmp = *this;
        --nIter;
        return tmp;
    }

    NeighborIteratorBase operator--(int) {
        --nIter;
        return *this;
    }

    bool operator==(const NeighborIteratorBase &rhs) const { return nIter == rhs.nIter; }

    bool operator!=(const NeighborIteratorBase &rhs) const { return !(nIter == rhs.nIter); }

    node operator*() const { return *nIter; }
};

/**
 * Class to iterate over the in/out neighbors of a node including the edge
 * weights. Values are std::pair<node, edgeweight>.
 */
class NeighborWeightIteratorBase {

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
    using self = NeighborWeightIteratorBase;

    NeighborWeightIteratorBase(std::vector<node>::const_iterator nodesIter,
                               std::vector<edgeweight>::const_iterator weightIter)
        : nIter(nodesIter), wIter(weightIter) {}

    /**
     * @brief WARNING: This contructor is required for Python and should not be used as the
     * iterator is not initialized.
     */
    NeighborWeightIteratorBase() {}

    NeighborWeightIteratorBase &operator++() {
        ++nIter;
        ++wIter;
        return *this;
    }

    NeighborWeightIteratorBase operator++(int) {
        const auto tmp = *this;
        ++(*this);
        return tmp;
    }

    NeighborWeightIteratorBase operator--() {
        --nIter;
        --wIter;
        return *this;
    }

    NeighborWeightIteratorBase operator--(int) {
        const auto tmp = *this;
        --(*this);
        return tmp;
    }

    bool operator==(const NeighborWeightIteratorBase &rhs) const {
        return nIter == rhs.nIter && wIter == rhs.wIter;
    }

    bool operator!=(const NeighborWeightIteratorBase &rhs) const { return !(*this == rhs); }

    std::pair<node, edgeweight> operator*() const { return std::make_pair(*nIter, *wIter); }
};

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_NEIGHBOR_ITERATORS_HPP_
