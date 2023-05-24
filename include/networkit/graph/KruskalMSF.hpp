/*
 * KruskalMSF.hpp
 *
 *  Created on: 03.09.2015
 *      Author: Henning
 */

#ifndef NETWORKIT_GRAPH_KRUSKAL_MSF_HPP_
#define NETWORKIT_GRAPH_KRUSKAL_MSF_HPP_

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/SpanningForest.hpp>

namespace NetworKit {

/**
 * Creates a minimum spanning tree for each connected component.
 * @ingroup graph
 */
class KruskalMSF final : public SpanningForest {
public:
    KruskalMSF(const Graph &G);

    /**
     * Computes for each component a minimum weight spanning tree
     * (or simply a spanning tree in unweighted graphs).
     * Uses Kruskal's algorithm.
     * Time complexity: sort(n) + n * inverse Ackermann(n, m).
     */
    void run() override;

    /**
     * @return Total edge-weight of the spanning forest.
     * Number of edges in forest if input graph is unweighted.
     */
    edgeweight getTotalWeight() const {
        assureFinished();
        if (G->isWeighted())
            return totalWeight;
        else
            return static_cast<edgeweight>(forest.numberOfEdges());
    }

private:
    edgeweight totalWeight = 0.0;
};

} /* namespace NetworKit */
#endif // NETWORKIT_GRAPH_KRUSKAL_MSF_HPP_
