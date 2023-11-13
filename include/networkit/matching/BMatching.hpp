#ifndef NETWORKIT_MATCHING_B_MATCHING_HPP_
#define NETWORKIT_MATCHING_B_MATCHING_HPP_

#include <set>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup b-matching
 */
class BMatching {

public:
    /**
     * Constructs a new BMatching.
     *
     * @param z Maximum number of nodes.
     * @param b b values
     */
    BMatching(const std::vector<count> &b, count z = 0);

    /**
     * Checks whether this is a proper b-matching.
     *
     * @param G
     * @return bool
     */
    bool isProper(const Graph &G) const;

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
     * @param G  A graph.
     * @return Number of edges in b-matching.
     */
    count size(const Graph &G) const;

    /**
     * Get total weight of edges in this b-matching.
     *
     * @param G
     * @return edgeweight
     */
    edgeweight weight(const Graph &G) const;

    const std::vector<std::set<node>> &getMatches() const;
    std::vector<count> getB() const;

protected:
    std::vector<std::set<node>> matches;
    const std::vector<count> b;
};

} /* namespace NetworKit */
#endif // NETWORKIT_MATCHING_B_MATCHING_HPP_
