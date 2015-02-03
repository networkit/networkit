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
 * Implements the Union Find data structure to maintain disjoint sets efficiently.
 */
class UnionFind {
private:
	std::vector<int> data;
public:
		
	/**
	 * Create a new set representation with not more the @max_element elements.
	 * Initially every element is in its own set.
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
	 * Find the representative to element @u
	 * @param u element
	 * @return representative of set containing @u
	 */	
	index find(index u);
	
	/**
	 *  Merge the two sets contain @u and @v
	 *  @param u element u
	 *  @param v element v
	 */
	void merge(index u, index v);

	/**
	 * Convert the Union Find data structure to a Partition
	 * @return Partition equivalent to the union find data structure
	 * */	
	Partition toPartition();
};
}
#endif
