#ifndef NETWORKIT_MATCHING_B_SUITOR_MATCHER_HPP_
#define NETWORKIT_MATCHING_B_SUITOR_MATCHER_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/matching/BMatcher.hpp>

namespace NetworKit {

/**
 * @ingroup matching
 * B-Suitor matching finding algorithm.
 */
class BSuitorMatcher final : public BMatcher {
public:
    /**
     * Computes a 1/2-approximate maximum weight b-matching of an undirected weighted Graph @c G
     * using the sequential b-Suitor algorithm published by Khan et al. in "Efficient Approximation
     * Algorithms For Weighted B-Matching", SIAM Journal on Scientific Computing, Vol. 38, Iss. 5
     * (2016).
     *
     * @param G An undirected graph.
     * @param b A value @a b that represents the max number of edges per vertex in the b-Matching.
     * Defaults to the ordinary 1-Matching.
     */
    BSuitorMatcher(const Graph &G, count b = 1); // TODO std::vector<count> &b

    ~BSuitorMatcher() override = default;

    /**
     * Runs the algorithm.
     */
    void run() override;

    /**
     * Iterates up to @a b times over the heaviest neighbors of node @a u and makes
     * them to suitors if eligible.
     *
     * @param u
     */
    void findSuitors(node u);

    /**
     * Finds the heaviest unmatched neighbor that @a y has not yet proposed to
     * if it exists. For equally weighted edges w(u, t), w(u, v) and t < v, w(u, t) is considered
     * smaller than w(u, v) to break ties.
     *
     * @param y
     * @return node
     */
    node findPreffered(node y);

    /**
     * Makes @a x a suitor of @a u and recursively calls itself for previous worse
     * suitors of @a u that got replaced with their new best match.
     *
     * @param u
     * @param x
     */
    void makeSuitor(node u, node x);

    /**
     * Inserts @a v into the sorted vector @a nodes[u] based on node weights. If there are
     * ties, the nodes are placed in lexicographic order. The size of @a nodes remains because the
     * last element gets removed.
     *
     * @param nodes
     * @param u
     * @param v
     */
    void sortInsert(std::vector<node> &nodes, node u, node v);

    /**
     * Removes @a u from @a nodes while maintaining its size. Subsequent nodes are shifted
     * and @c none is added to the end of @a nodes.
     *
     * @param nodes
     * @param u
     * @param v
     */
    void sortRemove(std::vector<node> &nodes, node u);

private:
    std::vector<std::vector<node>> suitors;
    std::vector<std::vector<node>> proposed;
    const int b;

    /**
     * Finds the index of the first @c none value in the list of @a nodes if present, otherwise @c
     * none.
     *
     * @param nodes
     * @param u
     * @return index
     */
    index findFirstFreeIndex(const std::vector<node> &nodes) const;

    /**
     * Finds the index of node @a x in the list of @a nodes if present, otherwise @c none.
     *
     * @param nodes
     * @param x
     * @return index
     */
    index findIndexOf(const std::vector<node> &nodes, node x) const;

    /**
     * Checks the symmetry of pairs of nodes. It must hold that v is in suitors(u) iff u is
     * in suitors(v).
     *
     */
    void checkSymmetry() const;
};
} // namespace NetworKit

#endif // NETWORKIT_MATCHING_B_SUITOR_MATCHER_HPP_
