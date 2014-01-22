/*
 * JaccardMeasure2.h
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef JACCARDMEASURE2_H_
#define JACCARDMEASURE2_H_

#include "DissimilarityMeasure2.h"

namespace NetworKit {

class JaccardMeasure2: public NetworKit::DissimilarityMeasure2 {

public:

	JaccardMeasure2();

	virtual ~JaccardMeasure2();

	virtual double getDissimilarity(Graph& G, Partition& first, Partition& second);
};

} /* namespace NetworKit */
#endif /* JACCARDMEASURE_H_ */
