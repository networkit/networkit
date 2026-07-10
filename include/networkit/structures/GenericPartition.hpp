/*
 * GenericPartition.hpp
 *
 *  Created on: 03.10.2013
 *      Author: cls
 */

#ifndef NETWORKIT_STRUCTURES_GENERIC_PARTITION_HPP_
#define NETWORKIT_STRUCTURES_GENERIC_PARTITION_HPP_

#include <cassert>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <networkit/Globals.hpp>

namespace NetworKit {

template <typename T>
concept IntegralValue = std::integral<T> && !std::same_as<std::remove_cvref_t<T>, bool>;

/**
 * @ingroup structures
 * Implements a partition of a set, i.e. a subdivision of the
 * set into disjoint subsets.
 */

template <IntegralValue IndexType>
class GenericPartition final {

    static constexpr IndexType noneIndex = std::numeric_limits<IndexType>::max();

public:
    GenericPartition();

    /**
     * Create a new partition data structure for @a z elements.
     *
     * @param[in] z maximum ValueType
     */
    GenericPartition(IndexType z);

    /**
     * Create a new partition data structure for @a z elements. Initialize each
     * entry to the default value. WARNING: this circumvents the standard
     * interface and may leave the object in an inconsistent state. Use only in
     * exceptional cases.
     *
     * @param[in] z maximum ValueType
     * @param[in] defaultValue
     */
    GenericPartition(IndexType z, IndexType defaultValue);

    GenericPartition(const std::vector<IndexType> &data);

    /* Updates the maximum index of the partition to @a newZ and sets all its
     * values to @a defaultValue.
     *
     * @param[in] newZ New maximum index of the partition. @param[in]
     * defaultValue Default value of all the elements in the partition.
     */
    void reset(IndexType newZ, IndexType defaultValue) {
        data.clear();
        data.resize(newZ, defaultValue);
        z = newZ;
        omega = 0;
    }

    /**
     *  Index operator.
     *
     *  @param[in] e an element
     */
    inline IndexType &operator[](IndexType e) { return this->data[e]; }

    /**
     * Index operator for const instances of this class.
     *
     * @param[in] e an element
     */
    inline const IndexType &operator[](IndexType e) const { return this->data[e]; }

    /**
     * Get the set (id) in which the element @a e is contained.
     *
     * @param e Index of element.
     * @return The index of the set in which @a e is contained.
     */
    inline IndexType subsetOf(IndexType e) const {
        assert(e < this->numberOfElements());
        return this->data[e];
    }

    /**
     * Extend the data structure and create a slot for one more element.
     * Initializes the entry to noneIndex and returns the index of the entry.
     */
    inline IndexType extend() {
        data.push_back(noneIndex);
        z++;
        assert(z == data.size()); //(data.size() - 1)
        return z - 1;
    }

    /**
     * Removes the entry for the given element
     * by setting it to noneIndex.
     */
    inline void remove(IndexType e) {
        assert(e < z);
        data[e] = noneIndex;
    }

    /**
     * Add a (previously unassigned) element @a e to the set @a s.
     *
     * @param s The index of the subset.
     * @param e The element to add.
     */
    inline void addToSubset(IndexType s, IndexType e) {
        assert(data[e] == noneIndex); // guarantee that element was unassigned
        assert(s <= omega);           // do not create new subset ids
        data[e] = s;
    }

    /**
     * Move the (previously assigned) element @a e to the set @a s.
     *
     * @param s The index of the subset.
     * @param e The element to move.
     */
    inline void moveToSubset(IndexType s, IndexType e) {
        assert(this->contains(e));
        assert(s <= omega); // do not create new subset ids
        data[e] = s;
    }

    /**
     * Creates a singleton set containing the element @a e.
     *
     * @param e The index of the element.
     */
    inline void toSingleton(IndexType e) { data[e] = newSubsetId(); }

    /**
     * Assigns every element to a singleton set.
     * Set id is equal to element id.
     */
    void allToSingletons();

    /**
     * Assigns every element to the same subset.
     * Set id is equal to zero.
     */
    void allToOnePartition();

    /**
     * Assigns the elements from both sets to a new set and returns the id of it.
     *
     * @param s Set to merge.
     * @param t Set to merge.
     * @return Id of newly created set.
     */
    IndexType mergeSubsets(IndexType s, IndexType t);

    /**
     * Sets an upper bound for the subset ids that CAN be assigned.
     *
     * @param[in] upper highest assigned subset ID + 1
     */
    inline void setUpperBound(IndexType upper) { this->omega = upper - 1; }

    /**
     * Return an upper bound for the subset ids that have been assigned.
     * (This is the maximum id + 1.)
     *
     * @return The upper bound.
     */
    inline IndexType upperBound() const { return omega + 1; }

    /**
     * Get a lower bound for the subset ids that have been assigned.
     *
     * @return The lower bound.
     */
    inline IndexType lowerBound() const { return 0; }

    /**
     * Change subset IDs to be consecutive, starting at 0.
     *
     * @param useTurbo Default: false. If set to true, a vector instead of a map to assign new ids
     * which results in a shorter running time but possibly a large space overhead.
     */
    void compact(bool useTurbo = false);

    /**
     * Check if partition assigns a valid subset to the element @a e.
     *
     * @param e The element.
     * @return @c true if the assigned subset is valid, @c false otherwise.
     */
    inline bool contains(IndexType e) const {
        // e is in the element index range and the entry is not empty
        return (e < z) && (data[e] != noneIndex);
    }

    /**
     * Check if two elements @a e1 and @a e2 belong to the same subset.
     *
     * @param e1 Element.
     * @param e2 Element.
     * @return @c true if @a e1 and @a e2 belong to same subset, @c false otherwise.
     */
    inline bool inSameSubset(IndexType e1, IndexType e2) const {
        assert(data[e1] != noneIndex);
        assert(data[e2] != noneIndex);
        return data[e1] == data[e2];
    }

    /**
     * Get a list of subset sizes. Indices do not necessarily correspond to subset ids.
     *
     * @return A vector of subset sizes.
     */
    std::vector<count> subsetSizes() const;

    /**
     * Get a map from subset id to size of the subset.
     *
     * @return A map from subset id to size of the subset.
     */
    std::map<IndexType, count> subsetSizeMap() const;

    /**
     * Get the members of the subset @a s.
     *
     * @param s The subset.
     * @return A set containing the members of @a s.
     */
    std::set<IndexType> getMembers(IndexType s) const;

    /**
     * @return number of elements in the partition.
     */
    inline IndexType numberOfElements() const {
        return z; // z is the maximum element id
    }

    /**
     * Get the current number of sets in this partition.
     *
     * @return The current number of sets.
     */
    count numberOfSubsets() const;

    /**
     * Get the actual vector representing the partition data structure.
     *
     * @return vector containing information about partitions.
     */
    const std::vector<IndexType> &getVector() const;

    /**
     * @return the subsets of the partition as a set of sets.
     */
    std::set<std::set<IndexType>> getSubsets() const;

    /**
     * Get the ids of nonempty subsets.
     *
     * @return A set of ids of nonempty subsets.
     */
    std::set<IndexType> getSubsetIds() const;

    /**
     * Set a human-readable identifier @a name for the instance.
     *
     * @param name The name.
     */
    inline void setName(std::string name) { this->name = std::move(name); }

    /**
     * Get the human-readable identifier.
     *
     * @return The name of this partition.
     */
    inline std::string getName() const { return this->name; }

    /**
     * Iterate over all entries (node, cluster id) and execute callback function @a func (lambda
     * closure).
     *
     * @param func Takes parameters <code>(node, index)</code>
     */
    template <typename Callback>
    void forEntries(Callback func) const;

    /**
     * Iterate over all entries (node, cluster id) in parallel and execute callback function @a
     * handle (lambda closure).
     *
     * @param handle Takes parameters <code>(node, index)</code>
     */
    template <typename Callback>
    void parallelForEntries(Callback handle) const;

private:
    IndexType z;     //!< maximum element index that can be mapped
    IndexType omega; //!< maximum subset index ever assigned
    std::vector<IndexType>
        data; //!< data container, indexed by element index, containing subset index
    std::string name;

    /**
     * Allocates and returns a new subset id.
     */
    inline IndexType newSubsetId() {
        IndexType s = ++omega;
        return s;
    }
};

template <IntegralValue IndexType>
template <typename Callback>
inline void GenericPartition<IndexType>::forEntries(Callback handle) const {
    for (IndexType e = 0; e < this->z; e++) {
        handle(e, data[e]);
    }
}

template <IntegralValue IndexType>
template <typename Callback>
inline void GenericPartition<IndexType>::parallelForEntries(Callback handle) const {
#pragma omp parallel for
    for (omp_index e = 0; e < static_cast<omp_index>(this->z); e++) {
        handle(static_cast<IndexType>(e), this->data[static_cast<IndexType>(e)]);
    }
}

using Partition = GenericPartition<index>;

} /* namespace NetworKit */

#include <networkit/structures/GenericPartitionImpl.hpp>

#endif // NETWORKIT_STRUCTURES_GENERIC_PARTITION_HPP_
