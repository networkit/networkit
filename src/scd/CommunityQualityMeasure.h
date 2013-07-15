/*
 * CommunityQualityMeasure.h
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef COMMUNITYQUALITYMEASURE_H_
#define COMMUNITYQUALITYMEASURE_H_

#include <unordered_set>
#include <math.h>
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * Abstract base class for all measures which evaluate the quality of a selectively found
 * community.
 */
class CommunityQualityMeasure {

protected:

	Graph* G;
	count degSum;

public:

	CommunityQualityMeasure(Graph& G);

	virtual ~CommunityQualityMeasure();

	virtual double getQuality(const std::unordered_set<node>& community) = 0;
};

class LocalModularity : public CommunityQualityMeasure{

public:

	LocalModularity(Graph&G);

	virtual ~LocalModularity();

	virtual double getQuality(const std::unordered_set<node>& community);
};

class Conduct : public CommunityQualityMeasure{

public:

	Conduct(Graph&G);

	virtual ~Conduct();

	virtual double getQuality(const std::unordered_set<node>& community);
};

} /* namespace NetworKit */
#endif /* COMMUNITYQUALITYMEASURE_H_ */
