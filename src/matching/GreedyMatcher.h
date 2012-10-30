/*
 * GreedyMatcher.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef GREEDYMATCHER_H_
#define GREEDYMATCHER_H_

#include "MaximalMatcher.h"

namespace EnsembleClustering {

class GreedyMatcher: public EnsembleClustering::MaximalMatcher {
public:
	GreedyMatcher();
	virtual ~GreedyMatcher();
};

} /* namespace EnsembleClustering */
#endif /* GREEDYMATCHER_H_ */
