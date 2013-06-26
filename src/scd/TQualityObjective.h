/*
 * TQualityObjective.h
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#ifndef TQUALITYOBJECTIVE_H_
#define TQUALITYOBJECTIVE_H_

#include "TGreedyCommunityExpansion.h"

namespace NetworKit {

class TQualityObjective {

public:

	TQualityObjective(const Graph& G, std::unordered_set<node>& community);

	virtual ~TQualityObjective();

	/**
	 * @param[in]	v	a candidate node
	 * @return the quality value achieved if we add the given node to the community
	 *
	 * Higher values are better.
	 */
	double getValue(node v);

public:

	const Graph* G;								//!< pointer to the graph
	std::unordered_set<node>* community;	//!< pointer to the current community

	count degSum; //!< degree sum of the graph needed
	count nBoundaryEdges; //!< current number of boundary edges
	count volume;	//!< current community volume
};

/**
 * LocalModularityM as a quality objective function
 */
class TLocalModularityM : public TQualityObjective {

public:

	TLocalModularityM(const Graph& G, std::unordered_set<node>& community);

	virtual ~TLocalModularityM();

	double getValue(node v);


};


/**
 * Conductance as a quality objective function. Unlike standard conductance,
 * higher values are better. This measure is defined as
 * $1 - conductance(C) = 1 - \frac{|B(C)|}{|\max \{vol (C), vol(\setminus \{ C \} )\}|}$
 */
class TConductance : public TQualityObjective {

public:

	TConductance(const Graph& G, std::unordered_set<node>& community);

	virtual ~TConductance();

	double getValue(node v);


};

} /* namespace NetworKit */
#endif /* TQUALITYOBJECTIVE_H_ */
