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

/**
 * @ingroup distance
 * Abstract base class for node distance measures.
 */
class NodeDistance {

protected:

	const Graph& G;

public:
	/**
	 * Constructs the NodeDistance class for the given Graph @a G.
	 *
	 * @param G The graph.
	 */
	NodeDistance(const Graph& G);

	/** Default destructor */
	virtual ~NodeDistance() = default;

	/**
	 * Perform preprocessing work. Needs to be called before distances are requested.
	 */
	virtual void preprocess() = 0;

	/**
	 * Return the distance between two nodes.
	 * The distance must be normed to return a distance between 0 and 1.
	 */
	virtual double distance(node u, node v) = 0;

	/**
	 * Returns the distances between all connected pairs of nodes.
	 * @return Vector containing the distances between all connected pairs of nodes.
	 */
	virtual std::vector<double> getEdgeScores() = 0;
};

} /* namespace NetworKit */
#endif /* NODEDISTANCE_H_ */
