/*
 * RandMeasure2.h
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NODESTRUCTURALRANDMEASURE2_H_
#define NODESTRUCTURALRANDMEASURE2_H_

#include "DissimilarityMeasure2.h"

namespace NetworKit {

/**
 * The node-structural Rand measure assigns a similarity value in [0,1]
 * to two partitions of a graph, by considering all pairs of nodes.
 */
class NodeStructuralRandMeasure2: public NetworKit::DissimilarityMeasure2 {

public:

	NodeStructuralRandMeasure2();

	virtual ~NodeStructuralRandMeasure2();

	virtual double getDissimilarity(Graph& G, Partition& first, Partition& second);

};

} /* namespace NetworKit */
#endif /* RANDMEASURE_H_ */
