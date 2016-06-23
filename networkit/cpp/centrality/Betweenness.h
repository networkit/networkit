/*
 * Betweenness.h
 *
 *  Created on: 19.02.2014
 *      Author: cls, ebergamini
 */

#ifndef BETWEENNESS_H_
#define BETWEENNESS_H_

#include "Centrality.h"

namespace NetworKit {

/**
 * @ingroup centrality
 */
class Betweenness: public NetworKit::Centrality {
public:
	/**
	 * Constructs the Betweenness class for the given Graph @a G. If the betweenness scores should be normalized,
	 * then set @a normalized to <code>true</code>. The run() method takes O(nm) time, where n is the number
	 * of nodes and m is the number of edges of the graph.
	 *
	 * @param G The graph.
	 * @param normalized Set this parameter to <code>true</code> if scores should be normalized in the interval [0,1].
	 * @param computeEdgeCentrality Set this parameter to <code>true</code> if edge betweenness should be computed as well.
	 */
	Betweenness(const Graph& G, bool normalized=false, bool computeEdgeCentrality=false);



	/**
	* Compute betweenness scores
	*
	*/
	void run() override;

	/*
	* Returns the maximum possible Betweenness a node can have in a graph with the same amount of nodes (=a star)
	*/
	double maximum() override;

};

} /* namespace NetworKit */

#endif /* BETWEENNESS_H_ */
