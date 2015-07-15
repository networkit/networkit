/*
 * ApproxCloseness.h
 *
 *  Created on: 16.06.2015
 *      Author: Arie Slobbe
 */

#ifndef APPROXCLOSENESS_H_
#define APPROXCLOSENESS_H_

#include "Centrality.h"


namespace NetworKit {

/**
 * @ingroup centrality
 * Approximation of closeness centrality according to algorithm described in
 * Eppstein, Wang: Fast Approximation of Centrality
 */
class ApproxCloseness: public NetworKit::Centrality {

public:

	/**
	 * The algorithm approximates the closeness of all nodes, by taking samples
	 * uniformly at random and solving the SSSP problem for each. More samples
	 * improves the accuracy of the approximation.
	 *
	 * @param	graph		input graph
	 * @param	nSamples	 user defined number of samples
	 * @param	normalized   normalize centrality values in interval [0,1] ?
	 */
	ApproxCloseness(const Graph& G, count nSamples, bool normalized=false);


	/**
	* Compute closeness scores parallel
	*
	*/
	void run() override;

	/*
	 * Returns the maximum possible Closeness a node can have in a graph with the same amount of nodes (=a star)
	 */
	double maximum();

private:

	count nSamples;

};

} /* namespace NetworKit */

#endif /* APPROXCLOSENESS_H_ */
