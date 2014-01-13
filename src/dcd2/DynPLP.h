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
	DynPLP(std::string prepStrategy="isolate", count theta=0);

	void update(std::vector<GraphEvent>& stream) override;

	Clustering detect() override;

private:
	std::string prepStrategy;

	count updateThreshold = 0;
	count nIterations = 0; //!< number of iterations in last run

	std::vector<bool> activeNodes;
};

} /* namespace NetworKit */

#endif /* DYNPLP_H_ */
