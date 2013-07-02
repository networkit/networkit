/*
 * NMIDistance.h
 *
 *  Created on: 30.04.2013
 *      Author: cls
 */

#ifndef NMIDISTANCE_H_
#define NMIDISTANCE_H_


#include "DissimilarityMeasure.h"
#include "DynamicNMIDistance.h"
#include "../overlap/HashingOverlapper.h"
#include "../overlap/RegionGrowingOverlapper.h"
#include "../auxiliary/MissingMath.h"
#include "../auxiliary/NumericTools.h"

namespace NetworKit {

/**
 * NMIDistance quantifies the dissimilarity between two clusterings using
 * Normalized Mutual Information.
 *
 */
class NMIDistance: public NetworKit::DissimilarityMeasure {

public:

	NMIDistance();

	virtual ~NMIDistance();

	virtual double getDissimilarity(Graph& G, Clustering& zeta, Clustering& eta);


};

} /* namespace NetworKit */
#endif /* NORMALIZEDMUTUALINFORMATION_H_ */
