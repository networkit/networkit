/*
 * DynCommunityDetector.h
 *
 *  Created on: 24.12.2013
 *      Author: cls
 */

#ifndef DYNCOMMUNITYDETECTOR_H_
#define DYNCOMMUNITYDETECTOR_H_

#include "../structures/Partition.h"
#include "../graph/Graph.h"


namespace NetworKit {

class DynCommunityDetector {

public:

	virtual void run(const Graph& G, Partition& communities) = 0;
};

} /* namespace NetworKit */

#endif /* DYNCOMMUNITYDETECTOR_H_ */
