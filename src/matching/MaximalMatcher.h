/*
 * MaximalMatcher.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef MAXIMALMATCHER_H_
#define MAXIMALMATCHER_H_

#include "Matcher.h"

namespace EnsembleClustering {

class MaximalMatcher: public EnsembleClustering::Matcher {
public:
	MaximalMatcher();
	virtual ~MaximalMatcher();
};

} /* namespace EnsembleClustering */
#endif /* MAXIMALMATCHER_H_ */
