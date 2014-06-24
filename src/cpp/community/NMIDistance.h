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
 * NMIDistance quantifies the dissimilarity between two clusterings using
 * Normalized Mutual Information.
 *
 */
class NMIDistance: public NetworKit::DissimilarityMeasure {

public:

	/** Default destructor */
	virtual ~NMIDistance();

	virtual double getDissimilarity(Graph& G, Partition& zeta, Partition& eta);

};

} /* namespace NetworKit */
#endif /* NORMALIZEDMUTUALINFORMATION_H_ */
