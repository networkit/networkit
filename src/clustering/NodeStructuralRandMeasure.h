/*
 * RandMeasure.h
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NODESTRUCTURALRANDMEASURE_H_
#define NODESTRUCTURALRANDMEASURE_H_

#include "DissimilarityMeasure.h"

namespace NetworKit {

class NodeStructuralRandMeasure: public NetworKit::DissimilarityMeasure {

public:

	NodeStructuralRandMeasure();

	virtual ~NodeStructuralRandMeasure();

	virtual double getDissimilarity(Graph& G, Clustering& first, Clustering& second);

};

} /* namespace NetworKit */
#endif /* RANDMEASURE_H_ */
