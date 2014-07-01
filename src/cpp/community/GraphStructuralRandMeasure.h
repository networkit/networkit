/*
 * GraphStructuralRandMeasure.h
 *
 *  Created on: 16.04.2013
 *      Author: cls
 */

#ifndef GRAPHSTRUCTURALRANDMEASURE_H_
#define GRAPHSTRUCTURALRANDMEASURE_H_

#include "DissimilarityMeasure.h"

namespace NetworKit {

/**
 * @ingroup community
 * The graph-structural Rand measure assigns a similarity value in [0,1]
 * to two partitions of a graph, by considering connected pairs of nodes.
 */
class GraphStructuralRandMeasure: public NetworKit::DissimilarityMeasure {

public:


	virtual double getDissimilarity(const Graph& G, const Partition& first, const Partition& second);

	//virtual double getDissimilarity(Graph& G, Clustering& zeta1, Graph& G2, Clustering& zeta2);
};

} /* namespace NetworKit */
#endif /* GRAPHSTRUCTURALRANDMEASURE_H_ */
