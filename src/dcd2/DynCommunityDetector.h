/*
 * DynCommunityDetector.h
 *
 *  Created on: 24.12.2013
 *      Author: cls
 */

#ifndef DYNCOMMUNITYDETECTOR_H_
#define DYNCOMMUNITYDETECTOR_H_

#include "../dynamics/GraphEvent.h"
#include "../clustering/Clustering.h"
#include "../graph/Graph.h"


namespace NetworKit {

class DynCommunityDetector {

public:

	DynCommunityDetector();

	virtual void attachGraph(Graph& G);

	virtual void update(std::vector<GraphEvent>& stream) = 0;

	virtual Clustering detect() = 0;

protected:

	Graph* G;
	Clustering zeta;
};

} /* namespace NetworKit */

#endif /* DYNCOMMUNITYDETECTOR_H_ */
