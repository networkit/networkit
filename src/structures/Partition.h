/*
 * Partition.h
 *
 *  Created on: 03.10.2013
 *      Author: cls
 */

#ifndef PARTITION_H_
#define PARTITION_H_

#include <cinttypes>
#include <set>
#include <vector>
#include <map>
#include <cassert>
#include <limits>



namespace NetworKit {

#define none std::numeric_limits<index>::max()	//!< absence of an entry

typedef uint64_t index;
typedef uint64_t count;

/**
 * Implements a partition of a set, i.e. a subdivision of the 
 * set into disjoint subsets.
 */
class Partition {

public:

	/**
	 * Create a new partition data structure for elements up to a maximum element index.
	 *
	 * @param[in]	z	maximum index
	 */
	Partition(index z);

	virtual ~Partition() = default;

	/**
	 *  Index operator.
	 *
	 *  @param[in]	e	an element
	 */
	inline index& operator [](const index& e) {
		return this->data[e];
	}
	/**
	 * Index operator for const instances of this class.
	 *
	 * @param[in]	e	an element
	 */
	inline const index& operator [](const index& e) const {
		return this->data[e];
	}

	/**
	 * Return the set (id) in which a element
	 * is contained.
	 */
	inline index subsetOf(index e) const {
		assert (e < this->numberOfElements());
		return this->data[e];
	}


	/**
	 * Add a (previously unassigned) element to a set
	 */
	void addToSubset(index s, index e);

	/**
	 * Move a (previously assigned) element to a set.
	 */
	void moveToSubset(index s, index e);

	/**
	 * Creates a singleton set containing the element.
	 */
	void toSingleton(index e);

	/**
	 * Assigns every element to a singleton set.
	 * Set id is equal to element id.
	 */
	void allToSingletons();

	/**
	 * Assigns the elements from both sets to a new set.
	 */
	void mergeSubsets(index s, index t);


	/**
	 * Check if partition is a 1-partition,
	 * i.e. every element is assigned to the same set.
	 */
	bool isOnePartition(const std::set<index>& elements);


	/**
	 * Check if partition is a singleton partition,
	 * i.e. every element is assigned to a different set.
	 */
	bool isSingletonPartition(const std::set<index>& elements) const;



	/**
	 * Return an upper bound for the subset ids that have been assigned.
	 *
	 * (This is the maximum id + 1.)
	 */
	index upperBound() const;

	/**
	 * Return a lower bound for the subset ids that have been assigned.
	 */
	index lowerBound() const;


	/**
	 * Change subset IDs to be consecutive, starting at 0.
	 */
	void compact();


	/**
	 * Check if partition assigns a valid subset to the element.
	 */
	bool contains(index e) const;


	/**
	 * Check if two elements belong to the same subset
	 */
	bool inSameSubset(index e1, index e2) const;


	/**
	 * Get a list of subset sizes. Indices do not necessarily correspond to subset ids.
	 */
	std::vector<count> subsetSizes() const;


	/**
	 * Get a map from subset id to size of the subset.
	 */
	std::map<index, count> subsetSizeMap() const;


	/**
	 * Get the members of a specific subset.
	 */
	std::set<index> getMembers(const index s) const;


	/**
	 * @return number of elements in the partition.
	 */
	count numberOfElements() const;


	/**
	 * Get the current number of sets in this partition.
	 */
	count numberOfSubsets() const;


	/**
	 * Iterate over all entries (node, cluster) and execute callback function (lambda closure).
	 */
	template<typename Callback> void forEntries(Callback func);


	/**
	 * Iterate over all entries (node, cluster) in parallel and execute callback function (lambda closure).
	 */
	template<typename Callback> void parallelForEntries(Callback handle);

	/**
	 * Iterate over all entries (node, cluster) and execute callback function (lambda closure).
	 */
	template<typename Callback> void forEntries(Callback func) const;


	/**
	 * Iterate over all entries (node, cluster) in parallel and execute callback function (lambda closure).
	 */
	template<typename Callback> void parallelForEntries(Callback handle) const;


private:
	index z;	//!< maximum element index that can be mapped
	index omega;	//!< maximum subset index ever assigned
	std::vector<index> data;  	//!< data container, indexed by element index, containing subset index

	/**
	 * Allocates and returns a new subset id.
	 */
	inline index newSubsetId() {
		omega++;
		index s = omega;
		return s;
	}
};

} /* namespace NetworKit */


template<typename Callback>
inline void NetworKit::Partition::forEntries(Callback handle) {
	for (index e = 0; e <= this->z; e += 1) {
		handle(e, data[e]);
	}

}

template<typename Callback>
inline void NetworKit::Partition::forEntries(Callback handle) const {
	for (index e = 0; e <= this->z; e += 1) {
		handle(e, data[e]);
	}
}

template<typename Callback>
inline void NetworKit::Partition::parallelForEntries(
		Callback handle) {
	#pragma omp parallel for
	for (index e = 0; e <= this->z; e += 1) {
		handle(e, data[e]);
	}
}


template<typename Callback>
inline void NetworKit::Partition::parallelForEntries(
		Callback handle) const {
	#pragma omp parallel for
	for (index e = 0; e <= this->z; e += 1) {
		handle(e, data[e]);
	}
}

#endif /* PARTITION_H_ */
