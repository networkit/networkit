/*
 * RandMeasure.h
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef RANDMEASURE_H_
#define RANDMEASURE_H_

#include "DissimilarityMeasure.h"

namespace NetworKit {

class RandMeasure: public NetworKit::DissimilarityMeasure {

public:

	RandMeasure();

	virtual ~RandMeasure();

	virtual double getDissimilarity(Graph& G, Clustering& first, Clustering& second);

};

} /* namespace NetworKit */
#endif /* RANDMEASURE_H_ */
