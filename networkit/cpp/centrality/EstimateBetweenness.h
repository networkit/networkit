/*
 *  EstimateBetweenness.h
 *
 *  Created on: 13.06.2014
 *      Author: Christian Staudt, Elisabetta Bergamini
 */

#ifndef ESTIMATEBETWEENNESS_H_
#define ESTIMATEBETWEENNESS_H_

#include "Centrality.h"


namespace NetworKit {

/**
 * @ingroup centrality
 * Estimation of betweenness centrality according to algorithm described in
 * Sanders, Geisberger, Schultes: Better Approximation of Betweenness Centrality
 */
class EstimateBetweenness: public NetworKit::Centrality {

public:

	/**
	 * The algorithm estimates the betweenness of all nodes, using weighting
	 * of the contributions to avoid biased estimation. The run() method takes O(m)
	 * time per sample, where  m is the number of edges of the graph.
	 *
	 * @param	graph		input graph
	 * @param	nSamples	 user defined number of samples
	 * @param	normalized   normalize centrality values in interval [0,1] ?
	 * @param	parallel_flag	if true, run in parallel with additional memory cost z + 3z * t
	 */
	 EstimateBetweenness(const Graph& G, count nSamples, bool normalized=false, bool parallel_flag=false);

	void run() override;


private:

	count nSamples;
	bool parallel_flag;

};

} /* namespace NetworKit */

#endif /* ESTIMATEBETWEENNESS_H_ */
