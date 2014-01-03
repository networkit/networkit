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

/*
 *
 */
class DynPLP: public NetworKit::DynCommunityDetector {
public:
	DynPLP(Graph& G, count theta);

	void process(std::vector<GraphEvent>& stream) override;

	Clustering retrieve() override;

private:

	count updateThreshold = 0;
	count nIterations = 0; //!< number of iterations in last run

	// TODO: std::vector<bool> activeNodes;
};

} /* namespace NetworKit */

#endif /* DYNPLP_H_ */
