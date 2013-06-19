/*
 * CommunityQualityMeasure.h
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef COMMUNITYQUALITYMEASURE_H_
#define COMMUNITYQUALITYMEASURE_H_

#include <unordered_set>

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Abstract base class for all measures which evaluate the quality of a selectively found
 * community.
 */
class CommunityQualityMeasure {

public:

	CommunityQualityMeasure();

	virtual ~CommunityQualityMeasure();

	virtual double getQuality(const std::unordered_set<node>& community, const Graph& G) = 0;
};

class LocalModularityL {

public:

	LocalModularityL();

	virtual ~LocalModularityL();

	virtual double getQuality(const std::unordered_set<node>& community, const Graph& G);
};

} /* namespace NetworKit */
#endif /* COMMUNITYQUALITYMEASURE_H_ */
