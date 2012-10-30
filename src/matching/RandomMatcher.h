/*
 * RandomMatcher.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef RANDOMMATCHER_H_
#define RANDOMMATCHER_H_

#include "MaximalMatcher.h"

namespace EnsembleClustering {

class RandomMatcher: public EnsembleClustering::MaximalMatcher {
public:
	RandomMatcher();
	virtual ~RandomMatcher();
};

} /* namespace EnsembleClustering */
#endif /* RANDOMMATCHER_H_ */
