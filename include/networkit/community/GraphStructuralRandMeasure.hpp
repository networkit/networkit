/*
 * GraphStructuralRandMeasure.h
 *
 *  Created on: 16.04.2013
 *      Author: cls
 */

#ifndef NETWORKIT_COMMUNITY_GRAPH_STRUCTURAL_RAND_MEASURE_HPP_
#define NETWORKIT_COMMUNITY_GRAPH_STRUCTURAL_RAND_MEASURE_HPP_

#include <networkit/community/DissimilarityMeasure.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * The graph-structural Rand measure assigns a similarity value in [0,1]
 * to two partitions of a graph, by considering connected pairs of nodes.
 */
class GraphStructuralRandMeasure: public DissimilarityMeasure {

public:


    virtual double getDissimilarity(const Graph& G, const Partition& first, const Partition& second);

    //virtual double getDissimilarity(Graph& G, Clustering& zeta1, Graph& G2, Clustering& zeta2);
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_GRAPH_STRUCTURAL_RAND_MEASURE_HPP_
