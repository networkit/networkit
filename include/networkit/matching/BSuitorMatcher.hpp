#ifndef NETWORKIT_MATCHING_B_SUITOR_MATCHER_HPP_
#define NETWORKIT_MATCHING_B_SUITOR_MATCHER_HPP_

#include <algorithm>
#include <networkit/graph/Graph.hpp>
#include <networkit/matching/BMatcher.hpp>

namespace NetworKit {

struct Node {
    node id;
    edgeweight weight;

    Node() = default;
    Node(node n, edgeweight w) : id(n), weight(w) {}

    bool operator==(const Node &other) const { return id == other.id && weight == other.weight; }
    bool operator!=(const Node &other) const { return id != other.id || weight != other.weight; }
};

struct Data {
    std::vector<Node> list;
    Node min; // (none, 0) if list still has free capacity
    edgeweight weight;
    count size;
    count max_size;

    Data() = default;

    Data(count b) {
        list.reserve(b);
        min = Node(none, 0);
        weight = 0;
        size = 0;
        max_size = b;
    }

    node popMinIfFull() {
        if (size < max_size) {
            return none;
        } else {
            auto ret = min.id;
            remove(min);
            return ret;
        }
    }

    void insert(const Node &v) {
        assert(size < max_size);

        list.emplace_back(v);
        size += 1;
        weight += v.weight;
        if (size == max_size && !list.empty()) {
            min = *std::min_element(list.begin(), list.end(), [](const Node &a, const Node &b) {
                if (a.weight == b.weight) {
                    return a.id > b.id;
                }
                return a.weight < b.weight;
            });
        }
    }

    void remove(const Node &u) {
        list.erase(
            std::remove_if(list.begin(), list.end(), [u](const Node &p) { return p.id == u.id; }),
            list.end());
        if (u.id == min.id) {
            min = Node(none, 0);
        }
        size -= 1;
        weight -= u.weight;
    }

    void sort() {
        std::sort(list.begin(), list.end(), [](const Node &u, const Node &v) {
            return (u.weight > v.weight || (u.weight == v.weight) && u.id < v.id);
        });
    }
};

struct Suitors : public Data {
    Suitors(count b) : Data(b) {}

    edgeweight getWeight() {
        edgeweight w = 0.0;

#pragma omp parallel for reduction(+ : w)
        for (auto e : list) {
            w += e.weight;
        }
        return w;
    }
};

struct Proposed : public Data {
    Proposed(count b) : Data(b) {}
};

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
     * @param b A vector of @a b values that represents the max number of edges per vertex @a v in
     * the b-Matching (b.at(v)).
     */
    BSuitorMatcher(const Graph &G, const std::vector<count> &b);

    /**
     * @param G An undirected graph.
     * @param b A value @a b that represents the max number of edges per vertex in the b-Matching.
     * Defaults to the ordinary 1-Matching.
     */
    BSuitorMatcher(const Graph &G, count b = 1);

    /**
     * @param G  An undirected graph.
     * @param path  A path to a file containing @a b values that represents the max number of edges
     * per vertex in the b-Matching.
     */
    BSuitorMatcher(const Graph &G, const std::string &path);

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
     * @return Node
     */
    Node findPreferred(node y);

    /**
     * Makes @a x a suitor of @a u and recursively calls itself for previous worse
     * suitors of @a u that got replaced with their new best match.
     *
     * @param u
     * @param w
     * @param x
     */
    void makeSuitor(node u, edgeweight w, node x);

    /**
     * Checks the symmetry of pairs of nodes. It must hold that v is in suitors(u) iff u is
     * in suitors(v).
     *
     */
    bool isSymmetrical() const;

    edgeweight getWeight() const;

private:
    std::vector<Suitors *> S;
    std::vector<Proposed *> T;
    const std::vector<count> b;

    /**
     * Reads values from a file at @a path into the vector of b-values.
     *
     * @param size
     * @param path
     * @return std::vector<count>
     */
    std::vector<count> readBValuesFromFile(count size, const std::string &path) const;
};
} // namespace NetworKit

#endif // NETWORKIT_MATCHING_B_SUITOR_MATCHER_HPP_
