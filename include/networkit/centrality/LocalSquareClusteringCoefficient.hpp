/*
 * LocalSquareClusteringCoefficient.hpp
 *
 *  Created on: 07.04.2022
 *      Author: tillahoffmann
 */

#ifndef NETWORKIT_CENTRALITY_LOCAL_SQUARE_CLUSTERING_COEFFICIENT_HPP_
#define NETWORKIT_CENTRALITY_LOCAL_SQUARE_CLUSTERING_COEFFICIENT_HPP_

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 */
class LocalSquareClusteringCoefficient : public Centrality {
public:
    /**
     * Constructs the LocalSquareClusteringCoefficient class for the given Graph @a G.
     *
     * @param G The graph.
     */
    LocalSquareClusteringCoefficient(const Graph &G);

    /**
     * Computes the local clustering coefficient on the graph passed in constructor.
     */
    void run() override;

    /**
     * Get the theoretical maximum of centrality score in the given graph.
     *
     * @return The maximum centrality score.
     */
    double maximum() override;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_LOCAL_SQUARE_CLUSTERING_COEFFICIENT_HPP_
