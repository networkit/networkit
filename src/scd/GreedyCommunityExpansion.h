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


/**
 * The Greedy Community Expansion algorithm.
 *
 * Greedily adds nodes from the shell to improve community quality.
 */
class GreedyCommunityExpansion: public NetworKit::SelectiveCommunityDetector {

public:

	class QualityObjective {
	public:
		Graph* G;								//!< pointer to the graph
		std::unordered_set<node>* community;	//!< pointer to the current community


	public:

		/**
		 * @param[in]	G	the graph
		 * @param[in]	community	the currently expanding community
		 */
		QualityObjective(Graph& G, std::unordered_set<node>& community);

		virtual ~QualityObjective();

		/**
		 * @param[in]	v	a candidate node
		 * @return the quality value achieved if we add the given node to the community
		 *
		 * Higher values are better.
		 */
		virtual double getValue(node v) = 0;
	};

	/**
	 * LocalModularityM as a quality objective function
	 */
	class LocalModularityM : public QualityObjective {

	public:

		LocalModularityM(Graph& G, std::unordered_set<node>& community);

		virtual ~LocalModularityM();

		virtual double getValue(node v);
	};

	/**
	 * Conductance as a quality objective function. Unlike standard conductance,
	 * higher values are better. This measure is defined as
	 * $1 - conductance(C) = 1 - \frac{|B(C)|}{|\max \{vol (C), vol(V�\setminus \{ C \} )\}|}$
	 */
	class Conductance : public QualityObjective {

	public:

		Conductance(Graph& G, std::unordered_set<node>& community);

		virtual ~Conductance();

		virtual double getValue(node v);
	};

	/**
	 * Acceptability measures quantify how likely a node from the community shell
	 * is to improve the community when it is included.
	 */
	class Acceptability {
	public:
		Graph* G;								//!< pointer to current graph
		std::unordered_set<node>* community;	//!< pointer to current community
		std::unordered_set<node>* shell;		//!< pointer to current shell
	public:

		/**
		 * @param[in]	G			pointer to current graph
		 * @param[in]	community	pointer to current community
		 * @param[in]	shell		pointer to current shell
		 */
		Acceptability(Graph& G, std::unordered_set<node>& community, std::unordered_set<node>& shell);

		virtual ~Acceptability();

		/**
		 * Get the acceptability value for a node. Higher values are better.
		 */
		virtual double getValue(node v) = 0;
	};


	class DummySimilarity: public Acceptability {

	public:

		DummySimilarity(Graph& G, std::unordered_set<node>& community, std::unordered_set<node>& shell);

		virtual ~DummySimilarity();

		virtual double getValue(node v);
	};

	/**
	 * Get the node cluster similarity value for a node.
	 *
	 * 	$\frac{|N(C) \cap N(v) |}{|N(C) \cup N(v)|}$
	 */
	class NodeClusterSimilarity: public Acceptability {

	public:

		NodeClusterSimilarity(Graph& G, std::unordered_set<node>& community, std::unordered_set<node>& shell);

		virtual ~NodeClusterSimilarity();

		virtual double getValue(node v);
	};




	/** Greedy Community Expansion **/

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



