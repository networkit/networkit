/*
 * RandomizedGreedyMatcher.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef RANDOMIZEDGREEDYMATCHER_H_
#define RANDOMIZEDGREEDYMATCHER_H_

#include "Matcher.h"

namespace EnsembleClustering {

class RandomizedGreedyMatcher: public EnsembleClustering::Matcher {
public:
	RandomizedGreedyMatcher();
	virtual ~RandomizedGreedyMatcher();
};

} /* namespace EnsembleClustering */
#endif /* RANDOMIZEDGREEDYMATCHER_H_ */
