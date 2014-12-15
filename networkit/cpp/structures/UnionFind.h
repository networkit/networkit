/*
 * UnionFind.h
 *
 *  Created on: 11.08.2014
 *      Author: Marcel Radermacher
 */

#ifndef UNIONFIND_H_
#define UNIONFIND_H_

#include <vector>
#include "../Globals.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 * @ingroup structures
 * Implements the Union Find Datastructure to maintain disjoint sets efficiently.
 */
class UnionFind {
private:
	std::vector<int> data;
public:
		
	/**
	 * Create a new set representation with not more the @max_element elements.
	 * Initialy every element is in its own set.
	 * @param max_element maximum number of elements 
	 */
	UnionFind(index max_element) : data(max_element) {
		allToSingletons();
	}

	/**
	 * Assigns every element to a singleton set.
	 * Set id is equal to element id.
	 */

	void allToSingletons();

	/**
	 * Find the represetive to element @u
	 * @param u element
	 * @return representive of set containing @u
	 */	
	index find(index u);
	
	/**
	 *  Merge the two sets contain @u and @v
	 *  @param u element u
	 *  @param v element v
	 */
	void merge(index u, index v);

	/**
	 * Convert the Union Find Datastructure to a Partition
	 * @return Partiion equivalent of the union find datastructure
	 * */	
	Partition toPartition();
};
}
#endif
