/*
 * BMatcher.hpp
 *
 *  Created on: 07.08.2023
 *      Author: Fabian Brandt-Tumescheit
 *              Frieda Gerharz
 */

#ifndef NETWORKIT_MATCHING_B_MATCHER_HPP_
#define NETWORKIT_MATCHING_B_MATCHER_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/matching/BMatching.hpp>

namespace NetworKit {

/**
 * @ingroup matching
 * Base class for b-matching algorithms.
 */
class BMatcher : public Algorithm {
protected:
    const Graph *G;
    BMatching bMatch;

public:
    /**
     * Constructs a new BMatcher.
     *
     * @param G The graph to compute the b-matching on.
     * @param b A vector of b-values for all nodes in the graph.
     */
    BMatcher(const Graph &G, const std::vector<count> &b);

    ~BMatcher() override = default;

    /**
     * Runs the b-matching algorithm on the stored graph.
     *
     */
    void run() override = 0;

    /**
     * Gets the b-matching for the stored graph.
     * Note: The b-matching is only valid, if the corresponding algorithm
     * has taken care of populating it properly.
     */
    BMatching getBMatching() const;
};

} /* namespace NetworKit */
#endif // NETWORKIT_MATCHING_B_MATCHER_HPP_
