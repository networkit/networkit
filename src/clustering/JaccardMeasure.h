/*
 * JaccardMeasure.h
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef JACCARDMEASURE_H_
#define JACCARDMEASURE_H_

#include "DissimilarityMeasure.h"

namespace NetworKit {

class JaccardMeasure: public NetworKit::DissimilarityMeasure {

public:

	JaccardMeasure();

	virtual ~JaccardMeasure();

	virtual double getDissimilarity(Graph& G, Clustering& first, Clustering& second);
};

} /* namespace NetworKit */
#endif /* JACCARDMEASURE_H_ */
