/*
 * ApproxBetweenness2.h
 *
 *  Created on: 13.06.2014
 *      Author: Christian Staudt, Elisabetta Bergamini
 */

#ifndef APPROXBETWEENNESS2_H_
#define APPROXBETWEENNESS2_H_

#include "Centrality.h"


namespace NetworKit {

/**
 * @ingroup centrality
 * Approximation of betweenness centrality according to algorithm described in
 * Sanders, Geisberger, Schultes: Better Approximation of Betweenness Centrality
 */
class ApproxBetweenness2: public NetworKit::Centrality {

public:

	/**
	 * The algorithm approximates the betweenness of all nodes, using weighting
	 * of the contributions to avoid biased estimation. The run() method takes O(m)
	 * time per sample, where  m is the number of edges of the graph.
	 *
	 * @param	graph		input graph
	 * @param	nSamples	 user defined number of samples
	 * @param	normalized   normalize centrality values in interval [0,1] ?
	 * @param	parallel_flag	if true, run in parallel with additional memory cost z + 3z * t
	 */
	ApproxBetweenness2(const Graph& G, count nSamples, bool normalized=false, bool parallel_flag=false);

	void run() override;


private:

	count nSamples;
	bool parallel_flag;

};

} /* namespace NetworKit */

#endif /* APPROXBETWEENNESS_H_ */
