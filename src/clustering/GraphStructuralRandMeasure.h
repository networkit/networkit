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

class GraphStructuralRandMeasure: public NetworKit::DissimilarityMeasure {

public:

	GraphStructuralRandMeasure();

	virtual ~GraphStructuralRandMeasure();

	virtual double getDissimilarity(Graph& G, Clustering& first, Clustering& second);
};

} /* namespace NetworKit */
#endif /* GRAPHSTRUCTURALRANDMEASURE_H_ */
