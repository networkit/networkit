/*
 * NodeDistance.h
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#ifndef NODEDISTANCE_H_
#define NODEDISTANCE_H_

#include "../graph/Graph.h"

namespace NetworKit {

class NodeDistance {

protected:

	const Graph& G;

public:

	NodeDistance(const Graph& G);

	virtual ~NodeDistance();

	/**
	 * Perform preprocessing work. Needs to be called before distances are requested.
	 */
	virtual void preprocess() = 0;

	/**
	 * Return the distance between two nodes.
	 */
	virtual void distance(node u, node v) = 0;
};

} /* namespace NetworKit */
#endif /* NODEDISTANCE_H_ */
