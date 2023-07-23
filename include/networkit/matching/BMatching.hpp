#ifndef NETWORKIT_MATCHING_B_MATCHING_HPP_
#define NETWORKIT_MATCHING_B_MATCHING_HPP_

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
     * @param b b value
     */
    BMatching(count z = 0, count b = 1);

    /**
     * Checks whether this is a proper b-matching.
     *
     * @param G
     * @return bool
     */
    bool isProper(const Graph &G) const;

    /**
     * Sets two nodes @a u and @a v as each others matching partners.
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
     * Get a vector of matched neighbors of @a v.
     * Elements with value @c none represent non existing nodes.
     *
     * @param v
     * @return std::vector<node>
     */
    std::vector<node> mates(node v) const;

    /**
     * Get total weight of edges in this b-matching.
     *
     * @param G
     * @return edgeweight
     */
    edgeweight weight(const Graph &G) const;

    std::vector<std::vector<node>> getMatrix() const;
    count getB() const;

protected:
    std::vector<std::vector<node>> data; //!< storage of b-matching nodes
    const count b;                       //!< b value

    /**
     * Finds the index of the first @c none value if present, otherwise @c none.
     *
     * @return index
     */
    index findFirstFreeIndex(node u) const;

    /**
     * Finds the index of node @a v in a list of matched nodes of @a u if present, otherwise @c
     * none.
     *
     * @param u
     * @param v
     * @return index
     */
    index findIndexOf(node u, node v) const;
};

} /* namespace NetworKit */
#endif // NETWORKIT_MATCHING_B_MATCHING_HPP_
