/*
 * Closeness.h
 *
 *  Created on: 03.10.2014
 *      Author: nemes
 */

#ifndef CLOSENESS_H_
#define CLOSENESS_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 */
class Closeness: public NetworKit::Centrality {
public:
	/**
	 * Constructs the Closeness class for the given Graph @a G. If the closeness scores should be normalized,
	 * then set @a normalized to <code>true</code>. The run() method takes O(nm) time, where n is the number
	 * of nodes and m is the number of edges of the graph. 
	 *
	 * @param G The graph.
	 * @param normalized Set this parameter to <code>true</code> if scores should be normalized in the interval [0,1].
	 * @param	checkConnectedness	turn this off if you know the graph is connected
	 *
	 * TODO: extend definition of closeness to disconnected graphs
	 */
	Closeness(const Graph& G, bool normalized=false, bool checkConnectedness=true);



	/**
	* Compute closeness scores parallel
	*
	*/
	void run() override;

	/*
	 * Returns the maximum possible Closeness a node can have in a graph with the same amount of nodes (=a star)
	 */
	double maximum() override;
};

} /* namespace NetworKit */

#endif /* CLOSENESS_H_ */
