/*
 * DynPLP.h
 *
 *  Created on: 03.01.2014
 *      Author: cls
 */

#ifndef DYNPLP_H_
#define DYNPLP_H_

#include "DynCommunityDetector.h"

namespace NetworKit {

typedef index label; // a label is the same as a community id


/*
 *
 */
class DynPLP: public NetworKit::DynCommunityDetector {
public:
	DynPLP(count theta=0);

	void process(std::vector<GraphEvent>& stream) override;

	Clustering retrieve() override;

private:

	count updateThreshold = 0;
	count nIterations = 0; //!< number of iterations in last run

	std::vector<bool> activeNodes;
};

} /* namespace NetworKit */

#endif /* DYNPLP_H_ */
