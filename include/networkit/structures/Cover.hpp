/*
 * Cover.hpp
 *
 *  Created on: 03.10.2013
 *      Author: cls
 */

#ifndef NETWORKIT_STRUCTURES_COVER_HPP_
#define NETWORKIT_STRUCTURES_COVER_HPP_

#include <map>
#include <set>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/structures/GenericPartition.hpp>

namespace NetworKit {

/**
 * @ingroup structures
 * Implements a cover of a set, i.e. an assignment of
 * its elements to possibly overlapping subsets.
 */
template <IntegralValue IndexType>
class GenericCover final {

public:
    /** Default constructor */
    GenericCover();

    /**
     * Create a new cover data structure for elements up to a maximum element index.
     *
     * @param[in] z maximum index
     */
    GenericCover(IndexType z);

    /**
     * Creates a new cover data structure which contains the given partition.
     *
     * @param[in] p The partition to construct the cover from
     */
    GenericCover(const GenericPartition<IndexType> &p);

    /**
     *  Index operator.
     *
     *  @param[in] e an element
     */
    inline std::set<IndexType> &operator[](const IndexType &e) { return this->data[e]; }
    /**
     * Index operator for const instances of this class.
     *
     * @param[in] e an element
     */
    inline const std::set<IndexType> &operator[](const IndexType &e) const { return this->data[e]; }

    /**
     * Return the ids of subsets in which the element @a e is contained.
     *
     * @param[in] e an element
     * @return A set of subset ids in which @a e is contained.
     */
    inline std::set<IndexType> subsetsOf(IndexType e) const {
        // TODO: assert (e < this->numberOfElements());
        return this->data[e];
    }

    /**
     * Check if cover assigns a valid subset to the element @a e.
     *
     * @param e an element.
     * @return @c true, if @a e is assigned to a valid subset, @c false otherwise.
     */
    bool contains(IndexType e) const;

    /**
     * Check if two elements @a e1 and @a e2 belong to the same subset.
     *
     * @param e1 an element.
     * @param e2 an element.
     * @return @c true, if @a e1 and @a e2 belong to the same subset, @c false otherwise.
     */
    bool inSameSubset(IndexType e1, IndexType e2) const;

    /**
     * Get the members of a specific subset @a s.
     *
     * @return The set of members of subset @a s.
     */
    std::set<IndexType> getMembers(IndexType s) const;

    /**
     * Add the (previously unassigned) element @a e to the set @a s.
     * @param[in] s a subset
     * @param[in] e an element
     */
    void addToSubset(IndexType s, IndexType e);

    /**
     * Remove the element @a e from the set @a s.
     * @param[in] s a subset
     * @param[in] e an element
     */
    void removeFromSubset(IndexType s, IndexType e);

    /**
     * Move the element @a e to subset @a s, i.e. remove it from all
     * other subsets and place it in the subset.
     *  @param[in] s a subset
     *  @param[in] e an element
     */
    void moveToSubset(IndexType s, IndexType e);

    /**
     * Creates a singleton set containing the element @a e and returns the index of the new set.
     * @param[in] e an element
     * @return The index of the new set.
     */
    IndexType toSingleton(IndexType e);

    /**
     * Assigns every element to a singleton set.
     * Set id is equal to element id.
     */
    void allToSingletons();

    /**
     * Assigns the elements from both sets to a new set.
     * @param[in] s a subset
     * @param[in] t a subset
     */
    void mergeSubsets(IndexType s, IndexType t);

    /**
     * Get an upper bound for the subset ids that have been assigned.
     * (This is the maximum id + 1.)
     *
     * @return An upper bound.
     */
    IndexType upperBound() const;

    /**
     * Get a lower bound for the subset ids that have been assigned.
     * @return A lower bound.
     */
    IndexType lowerBound() const;

    /**
     * Get a list of subset sizes. Indices do not necessarily correspond to subset ids.
     *
     * @return A list of subset sizes.
     */
    std::vector<count> subsetSizes() const;

    /**
     * Get a map from subset id to size of the subset.
     *
     * @return A map from subset id to size of the subset.
     */
    std::map<IndexType, count> subsetSizeMap() const;

    /**
     * Get the current number of sets in this cover.
     *
     * @return The number of sets in this cover.
     */
    count numberOfSubsets() const;

    /**
     * Get the current number of elements in this cover.
     *
     * @return The current number of elements.
     */
    count numberOfElements() const;

    /**
     * Add an additional element (node).
     */
    IndexType extend();

    /**
     * Get the ids of nonempty subsets.
     *
     * @return A set of ids of nonempty subsets.
     */
    std::set<IndexType> getSubsetIds() const;

    /**
     * Sets an upper bound for the subset ids that CAN be assigned.
     *
     * @param[in] upper highest assigned subset ID + 1
     */
    void setUpperBound(IndexType upper);

    /**
     * Iterate over all entries (node, subset ID of node) and execute callback function @a func
     * (lambda closure).
     *
     * @param func Takes parameters <code>(node, IndexType)</code>
     */
    template <typename Callback>
    void forEntries(Callback func) const;

    /**
     * Iterate over all entries (node, subset ID of node) in parallel and execute callback function
     * @a func (lambda closure).
     *
     * @param func Takes parameters <code>(node, IndexType)</code>
     */
    template <typename Callback>
    void parallelForEntries(Callback handle) const;

private:
    IndexType z;     //!< maximum element index that can be mapped
    IndexType omega; //!< maximum subset index ever assigned
                     //!< data container, indexed by element id, containing set of subset ids
    std::vector<std::set<IndexType>> data;

    /**
     * Allocates and returns a new subset id.
     */
    inline IndexType newSubsetId() {
        ++omega;
        IndexType s = omega;
        return s;
    }
};

template <IntegralValue IndexType>
template <typename Callback>
inline void GenericCover<IndexType>::forEntries(Callback handle) const {
    for (IndexType e = 0; e <= this->z; e += 1) {
        handle(e, data[e]);
    }
}

template <IntegralValue IndexType>
template <typename Callback>
inline void GenericCover<IndexType>::parallelForEntries(Callback handle) const {
#pragma omp parallel for
    for (omp_index e = 0; e <= static_cast<omp_index>(this->z); e += 1) {
        handle(static_cast<IndexType>(e), data[static_cast<IndexType>(e)]);
    }
}

using Cover = GenericCover<index>;

} /* namespace NetworKit */

#include <networkit/structures/CoverImpl.hpp>

#endif // NETWORKIT_STRUCTURES_COVER_HPP_
