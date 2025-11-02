#ifndef NETWORKIT_GRAPH_RANDOM_MAXIMUM_SPANNING_FOREST_HPP_
#define NETWORKIT_GRAPH_RANDOM_MAXIMUM_SPANNING_FOREST_HPP_

#include <limits>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/UnionFind.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/**
 * Computes a random maximum-weight spanning forest using Kruskal's algorithm by randomizing the
 * order of edges of the same weight.
 */
class RandomMaximumSpanningForest : public Algorithm {
public:
    /**
     * Initialize the random maximum-weight spanning forest algorithm, uses edge weights.
     *
     * @param G The input graph.
     */
    RandomMaximumSpanningForest(const Graph &G);

    /**
     * Initialize the random maximum-weight spanning forest algorithm using a vector of values, one
     * for each edge. These values are used as edge weights. Each value corresponds to the edge with
     * the given id.
     *
     * Note: Since the mapping relies on edge ids, this variant only works, if the graph has edge
     * ids. Call Graph::indexEdges() beforehand if necessary. The values are copied, the supplied
     * vector is not stored in the RandomMaximumSpanningForest object. If the graph is changed, i.e.
     * edges are removed, indexEdges should be called again to ensure that edge ids are contiguous.
     * Otherwise edgeValues entries for the missing edges can be arbitrary.
     *
     * @param G The input graph.
     * @param edgeValues The edge values to use, can be either of type edgeweight (double) or count
     * (uint64), internally all values are handled as double.
     */
    template <typename A>
    RandomMaximumSpanningForest(const Graph &G, const std::vector<A> &edgeValues);

    /**
     * Execute the algorithm. The algorithm is not parallel.
     */
    void run() override;

    /**
     * DEPRECATED: this function will no longer be supported in later releases. Use getIndicator()
     * instead.
     *
     * Get a boolean attribute that indicates for each edge if it is part of the calculated
     * maximum-weight spanning forest.
     *
     * This attribute is only calculated and can thus only be request if the supplied graph has edge
     * ids.
     *
     * @param move If the attribute shall be moved out of the algorithm instance.
     * @return The vector with the boolean attribute for each edge.
     */
    std::vector<bool> TLX_DEPRECATED(getAttribute(bool move = false));

    /**
     * Get a boolean indicator vector that indicates for each edge if it is part of any
     * maximum-weight spanning forest.
     *
     * This indicator vector is only calculated and can thus only be requested if the supplied graph
     * has edge ids.
     *
     * @param move If the indicator vector shall be moved out of the algorithm instance.
     * @return The vector with the boolean indicator for each edge.
     */
    std::vector<bool> getIndicator(bool move = false);

    /**
     * Checks if the edge (@a u, @a v) is part of the calculated maximum-weight spanning forest.
     *
     * @param u The first node of the edge to check
     * @param v The second node of the edge to check
     * @return If the edge is part of the calculated maximum-weight spanning forest.
     */
    bool inMSF(node u, node v) const;

    /**
     * Checks if the edge with the id @a eid is part of the calculated maximum-weight spanning
     * forest.
     *
     * @param eid The id of the edge to check.
     * @return If the edge is part of the calculated maximum-weight spanning forest.
     */
    bool inMSF(edgeid eid) const;

    /**
     * Gets the calculated maximum-weight spanning forest as graph.
     *
     * @param move If the graph shall be moved out of the algorithm instance.
     * @return The calculated maximum-weight spanning forest.
     */
    Graph getMSF(bool move = false);

private:
    struct weightedEdge {
        double attribute;
        node u;
        node v;
        edgeid eid;
        index rand;

        bool operator>(const weightedEdge &other) const {
            return (attribute > other.attribute)
                   || (attribute == other.attribute
                       && (rand > other.rand
                           || (rand == other.rand
                               && (u > other.u || (u == other.u && v > other.v)))));
        };
        weightedEdge(node u, node v, double attribute, edgeid eid = 0)
            : attribute(attribute), u(u), v(v), eid(eid), rand(Aux::Random::integer()){};
    };

    const Graph &G;
    std::vector<weightedEdge> weightedEdges;

    Graph msf;
    std::vector<bool> msfAttribute;

    bool hasWeightedEdges;
    bool hasMSF;
    bool hasAttribute;
};

template <typename A>
RandomMaximumSpanningForest::RandomMaximumSpanningForest(const Graph &G,
                                                         const std::vector<A> &edgeValues)
    : G(G), hasWeightedEdges(false), hasMSF(false), hasAttribute(false) {
    if (!G.hasEdgeIds()) {
        throw std::runtime_error("Error: Edges of G must be indexed for using edge values.");
    }

    weightedEdges.reserve(G.numberOfEdges());

    G.forEdges([&](node u, node v, edgeid eid) {
        weightedEdges.emplace_back(u, v, edgeValues[eid], eid);
    });

    INFO(weightedEdges.size(), " weighted edges saved");

    hasWeightedEdges = true;
}

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_RANDOM_MAXIMUM_SPANNING_FOREST_HPP_
