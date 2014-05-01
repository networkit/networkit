/*
 * DummySCD.h
 *
 *  Created on: 30.04.2013
 *      Author: cls
 */

#ifndef DUMMYSCD_H_
#define DUMMYSCD_H_

#include "SelectiveCommunityDetector.h"
#include "../graph/Graph.h"

namespace NetworKit {

class DummySCD : public SelectiveCommunityDetector {

public:

	DummySCD(const Graph& G);

	void run(std::set<unsigned int>& seeds) override;

};

} /* namespace NetworKit */
#endif /* DUMMYSCD_H_ */
