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


namespace NetworKit {

typedef uint64_t index;
typedef uint64_t count;

class Partition {

public:
	Partition(index z);
	virtual ~Partition();

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
	inline index setOf(index e) const {
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
	void toSingleton(index u);

	/**
	 * Assigns every element to a singleton set.
	 * Set id is equal to element id.
	 */
	void allToSingletons();

	/**
	 * Assigns the nodes from both sets to a new set.
	 */
	void mergeSubets(index s, index t);


	/**
	 * Check if partition is a 1-partition,
	 * i.e. every node is assigned to the same set.
	 */
	bool isOnePartition(const std::set<index>& elements);


	/**
	 * Check if partition is a singleton partition,
	 * i.e. every node is assigned to a different set.
	 */
	bool isSingletonPartition(const std::set<index>& elements) const;

	/**
	 * Get the current number of sets in this partition.
	 */
	count numberOfSubsets() const;


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
	 * Check if partition assigns a valid subset to the node.
	 */
	bool contains(index v) const;


	/**
	 * Check if two nodes belong to the same subset
	 */
	bool inSameSubset(index e1, index e2) const;


	bool equals(const Partition& other, const std::set<index>& elements) const;


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
	std::vector<index> getMembers(const index s) const;


	count numberOfElements() const;


private:
	std::vector<index> data;
	index nextset;	//!< next free set id for new set
	std::string name;
	index upperIdBound;	//!< upper bound for set ids


	/**
	 * Check if partition can hold a valid entry for the element because
	 * it is in the range mapped.
	 */
	bool isInRange(index v);

};

} /* namespace NetworKit */
#endif /* PARTITION_H_ */
