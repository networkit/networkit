/*
 * SelectiveSCAN.h
 *
 *  Created on: 14.06.2013
 *      Author: cls
 */

#ifndef SELECTIVESCAN_H_
#define SELECTIVESCAN_H_

#include "SelectiveCommunityDetector.h"
#include "../distmeasures/NodeDistance.h"

namespace NetworKit {

class SelectiveSCAN: public NetworKit::SelectiveCommunityDetector {



public:

	SelectiveSCAN(const Graph& G, NodeDistance& distMeasure, double epsilon=0.5, double mu=2);

	virtual ~SelectiveSCAN();

	virtual std::unordered_map<node, std::unordered_set<node> > run(std::unordered_set<node> seeds);


protected:


	void expandCore(node core, node label, std::unordered_set<node>* community,
					std::unordered_map<node, node>* nodesState, std::unordered_set<node>* candidates);

	std::pair<bool,std::unordered_set<node>> isCore(node u);


	double epsilon;
	double mu;

	NodeDistance* distMeasure;


};

} /* namespace NetworKit */
#endif /* SELECTIVESCAN_H_ */
