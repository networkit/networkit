/*
 * DummySCD.h
 *
 *  Created on: 30.04.2013
 *      Author: cls
 */

#ifndef DUMMYSCD_H_
#define DUMMYSCD_H_

#include "SelectiveCommunityDetector.h"

namespace NetworKit {

class DummySCD : public SelectiveCommunityDetector {

public:

	DummySCD(const Graph& G);

	void run(std::unordered_set<node> seeds) override;

	/** 
	 * @return a mapping from seed node to community (as a set of nodes)
	 */
	std::unordered_map<node, std::unordered_set<node> > getResult() override;

	/** 
	 * @return time in milliseconds spent on processing each seed node
	 */
	std::unordered_map<node, double> getTimings() override;

};

} /* namespace NetworKit */
#endif /* DUMMYSCD_H_ */
