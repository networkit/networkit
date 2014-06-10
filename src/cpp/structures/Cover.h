/*
 * Cover.h
 *
 *  Created on: 03.10.2013
 *      Author: cls
 */

#ifndef COVER_H_
#define COVER_H_

#include <cinttypes>
#include <set>
#include <vector>
#include <map>
#include <cassert>
#include <limits>


namespace NetworKit {

typedef uint64_t index;
typedef uint64_t count;

/**
 * Implements a cover of a set, i.e. an assignment of 
 * its elements to possibly overlapping subsets.
 */
class Cover {

public:
	Cover();

	/**
	 * Create a new cover data structure for elements up to a maximum element index.
	 *
	 * @param[in]	z	maximum index
	 */
	Cover(index z);

	virtual ~Cover() = default;


	/**
	 *  Index operator.
	 *
	 *  @param[in]	e	an element
	 */
	inline std::set<index>& operator [](const index& e) {
		return this->data[e];
	}
	/**
	 * Index operator for const instances of this class.
	 *
	 * @param[in]	e	an element
	 */
	inline const std::set<index>& operator [](const index& e) const {
		return this->data[e];
	}

	/**
	 * Return the ids of subsets in which a element
	 * is contained.
	 *
	 * 	 @param[in]	e	an element
	 */
	inline std::set<index> subsetsOf(index e) const {
		// TODO: assert (e < this->numberOfElements());
		return this->data[e];
	}



	/**
	 * Check if partition assigns a valid subset to the element.
	 */
	bool contains(index e) const;


	/**
	 * Check if two elements belong to the same subset
	 */
	bool inSameSubset(index e1, index e2) const;


	/**
	 * Get the members of a specific subset.
	 */
	std::set<index> getMembers(const index s) const;


	/**
	 * Add a (previously unassigned) element to a set
	 * 	@param[in]	s	a subset
	 *  @param[in]	e	an element
	 */
	void addToSubset(index s, index e);


	/**
	 * Move an element to a subset, i.e. remove it from all
	 * other subsets and place it in the subset.
	 * 	@param[in]	s	a subset
	 *  @param[in]	e	an element
	 */
	void moveToSubset(index s, index e);


	/**
	 * Creates a singleton set containing the element.
	 * @param[in]	e	an element
	 */
	index toSingleton(index e);


	/**
	 * Assigns every element to a singleton set.
	 * Set id is equal to element id.
	 */
	void allToSingletons();


	/**
	 * Assigns the elements from both sets to a new set.
	 * @param[in]	s	a subset
	 * @param[in]	t	a subset
	 */
	void mergeSubsets(index s, index t);


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
	 * Get a list of subset sizes. Indices do not necessarily correspond to subset ids.
	 */
	std::vector<count> subsetSizes() const;


	/**
	 * Get a map from subset id to size of the subset.
	 */
	std::map<index, count> subsetSizeMap() const;


	/**
	 * Get the current number of sets in this cover.
	 */
	count numberOfSubsets() const;

	/**
	 * Get the current number of elements in this cover.
	 */
	count numberOfElements() const;

	void setUpperBound(index upper);


	/**
	 * Iterate over all entries (node, set of subset IDs) and execute callback function (lambda closure).
	 */
	template<typename Callback> void forEntries(Callback func);


	/**
	 * Iterate over all entries (node, set of subset IDs) in parallel and execute callback function (lambda closure).
	 */
	template<typename Callback> void parallelForEntries(Callback handle);

	/**
	 * Iterate over all entries (node, set of subset IDs) and execute callback function (lambda closure).
	 */
	template<typename Callback> void forEntries(Callback func) const;


	/**
	 * Iterate over all entries (node, set of subset IDs) in parallel and execute callback function (lambda closure).
	 */
	template<typename Callback> void parallelForEntries(Callback handle) const;




private:

	index z;	//!< maximum element index that can be mapped
	index omega;	//!< maximum subset index ever assigned
	std::vector<std::set<index>> data;	//!< data container, indexed by element id, containing set of subset ids


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
inline void NetworKit::Cover::forEntries(Callback handle) {
	for (index e = 0; e <= this->z; e += 1) {
		handle(e, data[e]);
	}

}

template<typename Callback>
inline void NetworKit::Cover::forEntries(Callback handle) const {
	for (index e = 0; e <= this->z; e += 1) {
		handle(e, data[e]);
	}
}

template<typename Callback>
inline void NetworKit::Cover::parallelForEntries(
		Callback handle) {
	#pragma omp parallel for
	for (index e = 0; e <= this->z; e += 1) {
		handle(e, data[e]);
	}
}


template<typename Callback>
inline void NetworKit::Cover::parallelForEntries(
		Callback handle) const {
	#pragma omp parallel for
	for (index e = 0; e <= this->z; e += 1) {
		handle(e, data[e]);
	}
}


#endif /* COVER_H_ */
