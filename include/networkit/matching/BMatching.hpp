/*
 * BMatching.hpp
 *
 *  Created on: 07.08.2023
 *      Author: Fabian Brandt-Tumescheit
 *              Frieda Gerharz
 */

#ifndef NETWORKIT_MATCHING_B_MATCHING_HPP_
#define NETWORKIT_MATCHING_B_MATCHING_HPP_

#include <unordered_set>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup b-matching
 */
class BMatching {

public:
    /**
     * @brief WARNING: This contructor is required for Python and should not be used as the
     * BMatching is not correctly initialized.
     *
     */
    BMatching() = default;

    /**
     * Constructs a new BMatching.
     *
     * @param G The graph
     * @param b b values for all nodes
     */
    BMatching(const Graph &G, const std::vector<count> &b);

    /**
     * Checks whether this is a proper b-matching.
     *
     * @return bool
     */
    bool isProper() const;

    /**
     * Sets two nodes @a u and @a v as each others matching NodeMatches.
     *
     * @param u
     * @param v
     */
    void match(node u, node v);

    /**
     * Resets the two nodes @a u and @a v to unmatched.
     *
     * @param u
     * @param v
     */
    void unmatch(node u, node v);

    /**
     * Checks if node is unmatched.
     *
     * @param u
     * @return bool
     */
    bool isUnmatched(node u) const;

    /**
     * Checks if the two nodes @a u and @a v are matched together.
     *
     * @param u node.
     * @param v node.
     * @return bool
     */
    bool areMatched(node u, node v) const;

    /**
     * Get the number of edges in this b-matching.
     *
     * @return Number of edges in b-matching.
     */
    count size() const;

    /**
     * Get total weight of edges in this b-matching.
     *
     * @return edgeweight
     */
    edgeweight weight() const;

    /**
     * Retrieves a reference to the set of matches for each node.
     */
    const std::vector<std::unordered_set<node>> &getMatches() const;

    /**
     * Retrieves the b-value for each node.
     */
    const std::vector<count> &getB() const;

    /**
     * Removes all entries from the b-matching data structure
     */
    void reset();

protected:
    const Graph *G;
    std::vector<count> b;
    std::vector<std::unordered_set<node>> matches;
};

} /* namespace NetworKit */
#endif // NETWORKIT_MATCHING_B_MATCHING_HPP_
