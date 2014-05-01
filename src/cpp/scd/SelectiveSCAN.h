/*
 * SelectiveSCAN.h
 *
 *  Created on: 14.06.2013
 *      Author: cls, Yassien Marrakchi
 */

#ifndef SELECTIVESCAN_H_
#define SELECTIVESCAN_H_

#include "SelectiveCommunityDetector.h"
#include "../distmeasures/NodeDistance.h"

namespace NetworKit {
/**
 * Selective Clustering algorithm: selSCAN
 *
 */
class SelectiveSCAN: public NetworKit::SelectiveCommunityDetector {

public:

	SelectiveSCAN(const Graph& G, NodeDistance& distMeasure, double epsilon = 0.25, double mu = 3);

	virtual ~SelectiveSCAN();

	/**
	 * @param[in] seeds 	seed set
	 *
	 * Return the community and needed running time for each seed node
	 */
	void run(std::set<node>& seeds) override;


protected:

	/**
	 * Expand the community from a given core node
	 *
	 * @param[in] core			currently expanded core node
	 * @param[in] label			label of the current community
	 * @param[in] community		pointer to the currently expanded community
	 * @param[in] nodesState	pointer to the table of nodes assignment
	 * @param[in] candidates	pointer to the set of candidate nodes
	 */
	void expandCore(node core, node label, std::unordered_set<node>* community,
					std::unordered_map<node, node>* nodesState, std::unordered_set<node>* candidates);

	/**
	 * @param[in] node		considered node
	 *
	 * Verify if the considered node is a core node and return the set of close neighbors
	 */
	std::pair<bool,std::unordered_set<node>> isCore(node u);


	double epsilon;				//!< maximal distance for two close nodes
	double mu;					//!< minimal number of close neighbors to define a core
	NodeDistance* distMeasure; 	//!< pointer to the considered distance measure

};

} /* namespace NetworKit */
#endif /* SELECTIVESCAN_H_ */
