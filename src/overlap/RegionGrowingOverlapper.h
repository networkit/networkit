/*
 * RegionGrowingOverlapper.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef REGIONGROWINGOVERLAPPER_H_
#define REGIONGROWINGOVERLAPPER_H_

#include "Overlapper.h"

namespace EnsembleClustering {

class RegionGrowingOverlapper: public EnsembleClustering::Overlapper {
public:
	RegionGrowingOverlapper();
	virtual ~RegionGrowingOverlapper();
};

} /* namespace EnsembleClustering */
#endif /* REGIONGROWINGOVERLAPPER_H_ */
