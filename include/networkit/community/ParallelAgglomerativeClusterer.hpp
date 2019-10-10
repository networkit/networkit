/*
 * ParallelAgglomerativeClusterer.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu), Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#ifndef NETWORKIT_COMMUNITY_PARALLEL_AGGLOMERATIVE_CLUSTERER_HPP_
#define NETWORKIT_COMMUNITY_PARALLEL_AGGLOMERATIVE_CLUSTERER_HPP_

#include <networkit/community/CommunityDetectionAlgorithm.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * A parallel agglomerative community detection algorithm, maximizing modularity.
 */
class ParallelAgglomerativeClusterer: public CommunityDetectionAlgorithm {

public:
    /**
     * Constructor to the parallel agglomerative clusterer.
     *
     * @param[in]	G	input graph
     */
    ParallelAgglomerativeClusterer(const Graph& G);

    /**
     * Detect communities.
     */
    virtual void run();

    virtual std::string toString() const;
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_PARALLEL_AGGLOMERATIVE_CLUSTERER_HPP_
