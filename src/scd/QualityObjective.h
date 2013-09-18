/*
 * QualityObjective.h
 *
 *  Created on: 16.06.2013
 *      Author: cls, Yassine Marrakchi
 */

#ifndef QUALITYOBJECTIVE_H_
#define QUALITYOBJECTIVE_H_

#include "../graph/Graph.h"
#include <unordered_set>

namespace NetworKit {

/**
 * Quality objective quantifies the goodness of a given community.
 */
class QualityObjective {

public:

	/**
	 * @param[in]	G			the graph
	 * @param[in]	community	the currently expanding community
	 * @param[in]	boundary	the current boundary and the number of outgoing nodes pro node
	 */
	QualityObjective(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary);

	virtual ~QualityObjective();

	/**
	 * @param[in]	v	a candidate node
	 * @return the quality value achieved if we add the given node to the community
	 *
	 * Higher values are better.
	 */
	virtual std::vector<double> getValue(node v) = 0;

protected:
	const Graph* G;			//!< pointer to the graph
public:
	std::unordered_set<node>* community;		//!< pointer to the current community
	std::unordered_map<node,count>* boundary; 	//!< pointer to list of nodes laying  at the boundary and the
	count nBoundaryEdges; 	//!< current number of boundary edges
	count degSum; 			//!< degree sum of the graph needed
	count nNodes; 			//!< current number of nodes in the community
	count volume;			//!< current community volume
	count nInternEdges; 	//!< current number of intern edges
};


/**
 * Local modularity M as a quality objective function. Unlike standard local modularity M,
 * higher values are better. This measure is defined as the ratio of internal edges and
 * external edges
 */
class LocalModularityM : public QualityObjective {

public:

	LocalModularityM(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary);

	virtual ~LocalModularityM();

	virtual std::vector<double> getValue(node v);
};

/**
 * Local modularity L as a quality objective function.
 * Higher values are better.
 * This measure is defined as the ratio of internal density and external density.
 */
class LocalModularityL : public QualityObjective {

public:

	LocalModularityL(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary);

	virtual ~LocalModularityL();

	virtual std::vector<double> getValue(node v);
};
/**
 * Conductance as a quality objective function. Unlike standard conductance,
 * higher values are better. This measure is defined as
 * $1 - conductance(C) = 1 - \frac{|B(C)|}{|\max \{vol (C), vol(V�\setminus \{ C \} )\}|}$
 */
class ConductanceDistance : public QualityObjective {

public:

	ConductanceDistance(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary);

	virtual ~ConductanceDistance();

	virtual std::vector<double> getValue(node v);
};

} /* namespace NetworKit */
#endif /* QUALITYOBJECTIVE_H_ */
