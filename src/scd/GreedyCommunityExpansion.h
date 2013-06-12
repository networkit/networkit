/*
 * GreedyCommunityExpansion.h
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#ifndef GREEDYCOMMUNITYEXPANSION_H_
#define GREEDYCOMMUNITYEXPANSION_H_

#include <unordered_set>

#include "SelectiveCommunityDetector.h"


namespace NetworKit {

class GreedyCommunityExpansion: public NetworKit::SelectiveCommunityDetector {

public:

	class QualityObjective {

		/**
		 * @param[in]	G	the graph
		 * @param[in]	community	the currently expanding community
		 */
		QualityObjective(Graph& G, std::unordered_set<node>& community);

		virtual double getValue(node v) = 0;
	};


	class LocalModularity : public QualityObjective {

		// TODO: constructor

		virtual double getValue(node v);

	};


	class Acceptability {

		// TODO: constructor

		virtual double getValue(node v) = 0;

	};


	GreedyCommunityExpansion();

	virtual ~GreedyCommunityExpansion();

	/**
	 * @param[in]	G	input graph
	 * @param[in]	s	seed node
	 *
	 * @param[out]		community as a set of nodes
	 */
	virtual std::unordered_set<node> run(Graph& G, node s);
};

} /* namespace NetworKit */
#endif /* GREEDYCOMMUNITYEXPANSION_H_ */
