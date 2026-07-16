/*
 * GenericUnionFind.hpp
 *
 *  Created on: 11.08.2014
 *      Author: Marcel Radermacher
 *      Changed a bit by Henning Meyerhenke to reflect union by rank and path compression
 *        as taught in "Algorithms 1"
 */

#ifndef NETWORKIT_STRUCTURES_GENERIC_UNION_FIND_HPP_
#define NETWORKIT_STRUCTURES_GENERIC_UNION_FIND_HPP_

#include <cstddef>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/structures/GenericPartition.hpp>

namespace NetworKit {

/**
 * @ingroup structures
 * Implements the Union Find data structure to maintain disjoint sets efficiently.
 * Uses path compression and union by rank to achieve running time linear in
 * the number of elements times the inverse Ackermann function.
 */
template <IntegralValue IndexType>
class GenericUnionFind final {
    std::vector<IndexType> parent;
    std::vector<unsigned char> rank;

public:
    /**
     * Create a new set representation with not more than @p max_element elements.
     * Initially every element is in its own set.
     * @param max_element maximum number of elements
     */
    GenericUnionFind(IndexType max_element)
        : parent(static_cast<std::size_t>(max_element)),
          rank(static_cast<std::size_t>(max_element), 0) {
        allToSingletons();
    }

    /**
     * Assigns every element to a singleton set.
     * Set id is equal to element id.
     */
    void allToSingletons() {
        for (std::size_t i = 0; i < parent.size(); ++i) {
            parent[i] = static_cast<IndexType>(i);
        }
    }

    /**
     * Find the representative to element @u
     * @param u element
     * @return representative of set containing @u
     */
    IndexType find(IndexType u) {
        const auto uIndex = static_cast<std::size_t>(u);
        if (parent[uIndex] == u) {
            return u;
        }

        // recursion and path compression
        parent[uIndex] = find(parent[uIndex]);
        return parent[uIndex];
    }

    /**
     *  Merge the two sets contain @u and @v
     *  @param u element u
     *  @param v element v
     */
    void merge(IndexType u, IndexType v) {
        const auto setU = find(u);
        const auto setV = find(v);
        if (setU == setV) {
            return;
        }

        const auto setUIndex = static_cast<std::size_t>(setU);
        const auto setVIndex = static_cast<std::size_t>(setV);
        if (rank[setUIndex] < rank[setVIndex]) {
            parent[setUIndex] = setV;
        } else {
            parent[setVIndex] = setU;
            if (rank[setUIndex] == rank[setVIndex]) {
                ++rank[setUIndex];
            }
        }
    }

    /**
     * Convert the Union Find data structure to a Partition
     * @return Partition equivalent to the union find data structure
     * */
    GenericPartition<IndexType> toPartition() {
        GenericPartition<IndexType> partition(static_cast<IndexType>(parent.size()));
        partition.setUpperBound(static_cast<IndexType>(parent.size()));
        for (std::size_t e = 0; e < parent.size(); ++e) {
            const auto element = static_cast<IndexType>(e);
            partition.addToSubset(find(element), element);
        }
        return partition;
    }
};

using UnionFind = GenericUnionFind<index>;

} /* namespace NetworKit */

#endif // NETWORKIT_STRUCTURES_GENERIC_UNION_FIND_HPP_
