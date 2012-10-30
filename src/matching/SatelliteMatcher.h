/*
 * SatelliteMatcher.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef SATELLITEMATCHER_H_
#define SATELLITEMATCHER_H_

#include "RandomizedGreedyMatcher.h"

namespace EnsembleClustering {

class SatelliteMatcher: public EnsembleClustering::RandomizedGreedyMatcher {
public:
	SatelliteMatcher();
	virtual ~SatelliteMatcher();
};

} /* namespace EnsembleClustering */
#endif /* SATELLITEMATCHER_H_ */
