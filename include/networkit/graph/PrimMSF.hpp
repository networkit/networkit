
/*  PrimMST.hpp
 *
 *	Created on: 29.03.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#ifndef NETWORKIT_GRAPH_PRIM_MSF_HPP_
#define NETWORKIT_GRAPH_PRIM_MSF_HPP_
#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/SpanningForest.hpp>

namespace NetworKit {

/**
 * @class PrimMSF
 * @brief Computes a Minimum Spanning Forest using Prim's algorithm.
 *
 * For a given undirected graph, this algorithm grows a forest by greedily
 * adding the smallest edge between visited and unvisited nodes.
 * If the graph is connected, the result is a Minimum Spanning Tree (MST).
 * Otherwise, it yields a forest of MSTs for each connected component.
 */
class PrimMSF : public SpanningForest {
public:
    /**
     * @brief Initializes the PrimMSF algorithm for a given graph.
     *
     * The input graph must be undirected. It may be either weighted or unweighted.
     * In the case of an unweighted graph, each edge is treated as having unit weight.
     *
     * @param G The graph on which the minimum spanning forest will be computed.
     * @throws std::runtime_error if graph is not an undirected graph
     */
    PrimMSF(const Graph &G) : SpanningForest(G) {
        if (G.isDirected()) {
            throw std::runtime_error("The graph is not an undirected graph.");
        }
    }

    /**
     * @brief Runs Prim's algorithm to compute the minimum spanning forest.
     *
     * For each connected component of the graph, constructs a minimum spanning tree
     * and computes the total weight of the resulting forest.
     */
    void run() override;

    /**
     * @brief Returns the total weight of the minimum spanning forest.
     * @return The total weight of the minimum spanning forest.
     */
    edgeweight getTotalWeight() const {
        assureFinished();
        if (G->isWeighted())
            return totalWeight;
        return static_cast<edgeweight>(forest.numberOfEdges());
    }

private:
    static constexpr edgeweight infiniteWeight = std::numeric_limits<edgeweight>::max();
    edgeweight totalWeight = 0;
};
} // namespace NetworKit
#endif // NETWORKIT_GRAPH_PRIM_MSF_HPP_
