/*
 * APSP.hpp
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe
 */

#ifndef NETWORKIT_DISTANCE_APSP_HPP_
#define NETWORKIT_DISTANCE_APSP_HPP_

#include <memory>

#include <networkit/base/Algorithm.hpp>
#include <networkit/distance/SSSP.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Class for all-pair shortest path algorithm.
 */
class APSP : public Algorithm {

public:
    /**
     * Creates the APSP class for @a G.
     *
     * @param G The graph.
     */
    APSP(const Graph &G);

    ~APSP() override = default;

    /**
     * Computes the shortest paths from each node to all other nodes.
     * The algorithm is parallel.
     */
    void run() override;

    /**
     * Returns a vector of weighted distances between node pairs.
     *
     * @return The shortest-path distances from each node to any other node in
     * the graph.
     */
    const std::vector<std::vector<edgeweight>> &getDistances() const {
        assureFinished();
        return distances;
    }

    /**
     * Returns the distance from u to v or infinity if u and v are not
     * connected.
     *
     */
    edgeweight getDistance(node u, node v) const {
        assureFinished();
        return distances[u][v];
    }

protected:
    const Graph &G;
    std::vector<std::vector<edgeweight>> distances;
    std::vector<std::unique_ptr<SSSP>> sssps;
};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_APSP_HPP_
