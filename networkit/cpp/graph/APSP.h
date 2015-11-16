/*
 * APSP.h
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe
 */

#ifndef APSP_H_
#define APSP_H_

#include "Graph.h"
#include "../base/Algorithm.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Class for all-pair shortest path algorithm.
 */
class APSP: public Algorithm {

public:

	/**
	 * Creates the APSP class for @a G.
	 *
	 * @param G The graph.
	 */
	APSP(const Graph& G);

	virtual ~APSP() = default;

	/** Computes the shortest paths from each node to all other nodes. */
	void run() override;

	/**
	* @return string representation of algorithm and parameters.
	*/
	virtual std::string toString() const override;

	/**
	 * Returns a vector of weighted distances from the source node, i.e. the
 	 * length of the shortest path from the source node to any other node.
 	 *
 	 * @return The weighted distances from the source node to any other node in the graph.
	 */
	std::vector<std::vector<edgeweight> > getDistances() const { return distances;}


	/**
	 * Returns all shortest paths from source to @a t and an empty set if source and @a t are not connected.
	 *
	 */
	edgeweight getDistance(node u, node v) const { return distances[u][v];}

	/**
	* @return True if algorithm can run multi-threaded.
	*/
	virtual bool isParallel() const override { return true; }


protected:

	const Graph& G;
	std::vector<std::vector<edgeweight> > distances;
};

} /* namespace NetworKit */

#endif /* APSP_H_ */
