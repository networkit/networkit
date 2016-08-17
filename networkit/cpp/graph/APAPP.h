/*
 * APAPP.h
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe
 */

#ifndef APAPP_H_
#define APAPP_H_

#include "Graph.h"
#include "../base/Algorithm.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Class for all-pair algebraic path algorithm.
 */
template <class T>
class APAPP: public Algorithm {

public:

	/**
	 * Creates the APAPP class for @a G.
	 *
	 * @param G The graph.
	 */
	APAPP(const Graph& G, std::vector<T> edgeWeights={}) : Algorithm(), G(G), edgeWeights(edgeWeights) {
    }

	virtual ~APAPP() = default;

	/** Computes the shortest paths from each node to all other nodes. */
	virtual void run() override = 0;

	/**
	* @return string representation of algorithm and parameters.
	*/
	virtual std::string toString() const override {
        return "All-pairs Algebraic Path Algorithm";
    }

	/**
	 * Returns a vector of weighted distances from the source node, i.e. the
 	 * length of the shortest path from the source node to any other node.
 	 *
 	 * @return The weighted distances from the source node to any other node in the graph.
	 */
	std::vector<std::vector<T> > getDistances() const { return distances;}


	/**
	 * Returns all shortest paths from source to @a t and an empty set if source and @a t are not connected.
	 *
	 */
	T getDistance(node u, node v) const { return distances[u][v];}

	/**
	* @return True if algorithm can run multi-threaded.
	*/
	virtual bool isParallel() const override { return true; }


protected:

	const Graph& G;
	std::vector<std::vector<T> > distances;
    std::vector<T> edgeWeights;
};

} /* namespace NetworKit */

#endif /* APAPP_H_ */
