/*
 *
 * CommunityQualityMeasure.h
 *
 * Created on: 18.06.2013
 * Author: cls, Yassine Marrakchi
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

	Graph* G; //!< pointer to graph
	count degSum; //!< graphs's volume

public:

	/**
	 * @param[in] G		pointer to current graph
	 */
	CommunityQualityMeasure(Graph& G);

	virtual ~CommunityQualityMeasure();

	/**
	 * Get the acceptability value for a node. Higher values are better.
	 *
	 * @param[in] community 	pointer to the evaluated community
	 *
	 * return the objective value of the considered community
	 */
	virtual double getQuality(const std::unordered_set<node>& community) = 0;
};

/**
 * Implementation of local modularity L defined as the internal density
 * divided by the external density.
 * Higher values are better
 *
 */
class LocalModularity : public CommunityQualityMeasure{

public:

	LocalModularity(Graph&G);

	virtual ~LocalModularity();

	virtual double getQuality(const std::unordered_set<node>& community);
};

/**
 * Implementation of conductance defined as the ratio of boundary size and
 * the minimum of vol(C) and vol(V/C)
 *
 * Lower Values are better
 */
class Conduct : public CommunityQualityMeasure{

public:

	Conduct(Graph&G);

	virtual ~Conduct();

	virtual double getQuality(const std::unordered_set<node>& community);
};

} /* namespace NetworKit */
#endif /* COMMUNITYQUALITYMEASURE_H_ */
