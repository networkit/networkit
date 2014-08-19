/*
 * ClusteredRandomGraphGenerator.h
 *
 *  Created on: 28.02.2014
 *      Author: cls
 */

#ifndef CLUSTEREDRANDOMGRAPHGENERATOR_H_
#define CLUSTEREDRANDOMGRAPHGENERATOR_H_

#include "StaticGraphGenerator.h"

namespace NetworKit {

class ClusteredRandomGraphGenerator: public NetworKit::StaticGraphGenerator {
public:
	/**
	 * Creates a clustered random graph:
	 *
	 * @param[in]	n	number of nodes
	 * @param[in]	k	number of clusters
	 * @param[in]	pin		intra-cluster edge probability
	 * @param[in]	pout	inter-cluster edge probability
	 */
	ClusteredRandomGraphGenerator(count n, count k, double pin, double pout);

	Graph generate() override;

private:

	count n;
	count k;
	double pin;
	double pout;
};

} /* namespace NetworKit */

#endif /* CLUSTEREDRANDOMGRAPHGENERATOR_H_ */
