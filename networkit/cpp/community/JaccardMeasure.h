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

/**
 * @ingroup community
 */
class JaccardMeasure: public NetworKit::DissimilarityMeasure {

public:

	double getDissimilarity(const NetworKit::Graph &G, const NetworKit::Partition &zeta, const NetworKit::Partition &eta) override;
};

} /* namespace NetworKit */
#endif /* JACCARDMEASURE_H_ */
