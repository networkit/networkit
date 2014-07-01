/*
 * NMIDistance.h
 *
 *  Created on: 30.04.2013
 *      Author: cls
 */

#ifndef NMIDISTANCE_H_
#define NMIDISTANCE_H_


#include "DissimilarityMeasure.h"

namespace NetworKit {

/**
 * @ingroup community
 * NMIDistance quantifies the dissimilarity between two clusterings using
 * Normalized Mutual Information.
 *
 */
class NMIDistance: public NetworKit::DissimilarityMeasure {

public:


	virtual double getDissimilarity(const Graph& G, const Partition& zeta, const Partition& eta);

};

} /* namespace NetworKit */
#endif /* NORMALIZEDMUTUALINFORMATION_H_ */
