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

	TQualityObjective(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary);

	virtual ~TQualityObjective();

	/**
	 * @param[in]	v	a candidate node
	 * @return the quality value achieved if we add the given node to the community
	 *
	 * Higher values are better.
	 */
	std::vector<double> getValue(node v);

protected:
	const Graph* G;								//!< pointer to the graph
public:
	std::unordered_set<node>* community;	//!< pointer to the current community
	std::unordered_map<node,count>* boundary;
	count nBoundaryEdges; //!< current number of boundary edges
	count degSum; //!< degree sum of the graph needed
	count nNodes; //!< current number of nodes in the community
	count volume;	//!< current community volume
	count nInternEdges; //!< current number of intern edges
};

/**
 * LocalModularityM as a quality objective function
 */
class TLocalModularityM : public TQualityObjective {

public:

	TLocalModularityM(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary);

	virtual ~TLocalModularityM();

	std::vector<double> getValue(node v);


};

class TLocalModularityL : public TQualityObjective {

public:

	TLocalModularityL(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary);

	virtual ~TLocalModularityL();

	std::vector<double> getValue(node v);
};

/**
 * Conductance as a quality objective function. Unlike standard conductance,
 * higher values are better. This measure is defined as
 * $1 - conductance(C) = 1 - \frac{|B(C)|}{|\max \{vol (C), vol(\setminus \{ C \} )\}|}$
 */
class TConductance : public TQualityObjective {

public:

	TConductance(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary);

	virtual ~TConductance();

	std::vector<double> getValue(node v);


};

} /* namespace NetworKit */
#endif /* TQUALITYOBJECTIVE_H_ */
