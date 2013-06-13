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
	public:
		Graph G;
		std::unordered_set<node>* community;

		/**
		 * @param[in]	G	the graph
		 * @param[in]	community	the currently expanding community
		 */
	public:
		QualityObjective(Graph& G, std::unordered_set<node>* community);
		virtual ~QualityObjective();
		virtual double getValue(node v) = 0;
	};

	/**
	 * LocalModularityM as a quality objective function
	 */
	class LocalModularityM : public QualityObjective {
	public:
		LocalModularityM(Graph& G, std::unordered_set<node>* community);
		virtual ~LocalModularityM();
		virtual double getValue(node v);
	};

	class Conductance : public QualityObjective {
		public:
		Conductance(Graph& G, std::unordered_set<node>* community);
		virtual ~Conductance();
		virtual double getValue(node v);
		};

	/**
	 * Acceptability measures quantify how likely a node from the community shell
	 * is to improve the community when it is included.
	 */
	class Acceptability {
	public:
		Graph G;
		std::unordered_set<node>* community;
		std::unordered_set<node>* shell;
	public:
		Acceptability(Graph& G, std::unordered_set<node>* community, std::unordered_set<node>* shell);
		virtual ~Acceptability();
		virtual double getValue(node v) = 0;
	};


	class NodeClusterSimilarity: public Acceptability {
	public:
		NodeClusterSimilarity(Graph& G, std::unordered_set<node>* community, std::unordered_set<node>* shell);
		virtual ~NodeClusterSimilarity();
		virtual double getValue(node v);
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



