/*
 * NeighborIterators.hpp
 *
 *  Created on: 16.09.2024
 */

#ifndef NETWORKIT_GRAPH_NEIGHBOR_ITERATORS_HPP_
#define NETWORKIT_GRAPH_NEIGHBOR_ITERATORS_HPP_

namespace NetworKit {

/**
 * Class to iterate over the in/out neighbors of a node.
 */
template <typename containerType>
class NeighborIteratorBase {

    typename containerType::const_iterator elementIter;

public:
    // The value type of the neighbors (i.e. nodes). Returned by
    // operator*().
    using value_type = typename containerType::value_type;

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

    NeighborIteratorBase(typename containerType::const_iterator elementIter)
        : elementIter(elementIter) {}

    /**
     * @brief WARNING: This contructor is required for Python and should not be used as the
     * iterator is not initialized.
     */
    NeighborIteratorBase() {}

    NeighborIteratorBase &operator++() {
        ++elementIter;
        return *this;
    }

    NeighborIteratorBase operator++(int) {
        const auto tmp = *this;
        ++elementIter;
        return tmp;
    }

    NeighborIteratorBase operator--() {
        const auto tmp = *this;
        --elementIter;
        return tmp;
    }

    NeighborIteratorBase operator--(int) {
        --elementIter;
        return *this;
    }

    bool operator==(const NeighborIteratorBase &rhs) const {
        return elementIter == rhs.elementIter;
    }

    bool operator!=(const NeighborIteratorBase &rhs) const {
        return !(elementIter == rhs.elementIter);
    }

    value_type operator*() const { return *elementIter; }
};

/**
 * Class to iterate over the in/out neighbors of a containerType, including the
 * weights given by weightContainerType.
 */
template <typename containerType, typename weightContainerType>
class NeighborWeightIteratorBase {

    typename containerType::const_iterator elementIter;
    typename weightContainerType::const_iterator weightIter;

public:
    // The value type of the neighbors (i.e. nodes). Returned by
    // operator*().
    using value_type =
        std::pair<typename containerType::value_type, typename weightContainerType::value_type>;

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

    NeighborWeightIteratorBase(typename containerType::const_iterator elementIter,
                               typename weightContainerType::const_iterator weightIter)
        : elementIter(elementIter), weightIter(weightIter) {}

    /**
     * @brief WARNING: This contructor is required for Python and should not be used as the
     * iterator is not initialized.
     */
    NeighborWeightIteratorBase() {}

    NeighborWeightIteratorBase &operator++() {
        ++elementIter;
        ++weightIter;
        return *this;
    }

    NeighborWeightIteratorBase operator++(int) {
        const auto tmp = *this;
        ++(*this);
        return tmp;
    }

    NeighborWeightIteratorBase operator--() {
        --elementIter;
        --weightIter;
        return *this;
    }

    NeighborWeightIteratorBase operator--(int) {
        const auto tmp = *this;
        --(*this);
        return tmp;
    }

    bool operator==(const NeighborWeightIteratorBase &rhs) const {
        return elementIter == rhs.elementIter && weightIter == rhs.weightIter;
    }

    bool operator!=(const NeighborWeightIteratorBase &rhs) const { return !(*this == rhs); }

    std::pair<typename containerType::value_type, typename weightContainerType::value_type>
    operator*() const {
        return std::make_pair(*elementIter, *weightIter);
    }
};

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_NEIGHBOR_ITERATORS_HPP_
