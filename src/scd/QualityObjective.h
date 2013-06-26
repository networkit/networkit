/*
 * QualityObjective.h
 *
 *  Created on: 16.06.2013
 *      Author: Yassine Marrakchi
 */

#ifndef QUALITYOBJECTIVE_H_
#define QUALITYOBJECTIVE_H_

#include "../graph/Graph.h"
#include <unordered_set>

namespace NetworKit {

class QualityObjective {

public:

		/**
		 * @param[in]	G	the graph
		 * @param[in]	community	the currently expanding community
		 */
	QualityObjective(const Graph& G, std::unordered_set<node>& community);

	virtual ~QualityObjective();

		/**
		 * @param[in]	v	a candidate node
		 * @return the quality value achieved if we add the given node to the community
		 *
		 * Higher values are better.
		 */
	virtual double getValue(node v) = 0;

protected:
	const Graph* G;								//!< pointer to the graph
	std::unordered_set<node>* community;	//!< pointer to the current community
};


/**
 * LocalModularityM as a quality objective function
 */
class LocalModularityM : public QualityObjective {

public:

	LocalModularityM(const Graph& G, std::unordered_set<node>& community);

	virtual ~LocalModularityM();

	virtual double getValue(node v);
};

class LocalModularityL : public QualityObjective {

public:

	LocalModularityL(const Graph& G, std::unordered_set<node>& community);

	virtual ~LocalModularityL();

	virtual double getValue(node v);
};
/**
 * Conductance as a quality objective function. Unlike standard conductance,
 * higher values are better. This measure is defined as
 * $1 - conductance(C) = 1 - \frac{|B(C)|}{|\max \{vol (C), vol(V�\setminus \{ C \} )\}|}$
 */
class Conductance : public QualityObjective {

public:

	Conductance(const Graph& G, std::unordered_set<node>& community);

	virtual ~Conductance();

	virtual double getValue(node v);

public: // TODO: make this protected

	count degSum; //!< degree sum of the graph needed
	count nBoundaryEdges; //!< current number of boundary edges
	count volume;	//!< current community volume

};

} /* namespace NetworKit */
#endif /* QUALITYOBJECTIVE_H_ */
