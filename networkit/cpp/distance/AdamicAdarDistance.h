/*
 * AdamicAdarDistance.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#ifndef ADAMICADARDISTANCE_H_
#define ADAMICADARDISTANCE_H_

#include "../graph/Graph.h"
#include "NodeDistance.h"

namespace NetworKit {

/**
 * @ingroup distance
 * An implementation of the Adamic Adar distance measure.
 */
class AdamicAdarDistance : public NodeDistance {

protected:
	std::vector<double> aaDistance; //result vector

	void removeNode(Graph& graph, node u);

public:

	/**
	 * @param G The graph.
	 */
	AdamicAdarDistance(const Graph& G);

	/**
	 * Computes the Adamic Adar distances of all connected pairs of nodes.
	 * REQ: Needs to be called before distance() and getEdgeScores() deliver meaningful results!
	 */
	virtual void preprocess();

	/**
	 * Returns the Adamic Adar distance between node @a u and node @a v.
	 * @return Adamic Adar distance between the two nodes.
	 */
	virtual double distance(node u, node v);

	/**
	 * Returns the Adamic Adar distances between all connected nodes.
	 * @return Vector containing the Adamic Adar distances between all connected pairs of nodes.
	 */
	virtual std::vector<double> getEdgeScores();

};

} /* namespace NetworKit */

#endif /* ADAMICADARDISTANCE_H_ */
