/*
 * ApproxBetweenness.h
 *
 *  Created on: 09.04.2014
 *      Author: cls
 */

#ifndef APPROXBETWEENNESS_H_
#define APPROXBETWEENNESS_H_

#include "Centrality.h"

namespace NetworKit {


/**
 * @ingroup centrality
 * Approximation of betweenness centrality according to algorithm described in
 * Matteo Riondato and Evgenios M. Kornaropoulos: Fast Approximation of Betweenness Centrality through Sampling
 */
class ApproxBetweenness: public NetworKit::Centrality {

public:

	/**
	 * The algorithm approximates the betweenness of all vertices so that the scores are
	 * within an additive error @a epsilon with probability at least (1- @a delta).
	 * The values are normalized by default.
	 *
	 * @param	G			the graph
	 * @param	epsilon		maximum additive error
	 * @param	delta		probability that the values are within the error guarantee
	 * @param	diameterSamples		if 0, use the possibly slow estimation of the vertex diameter which definitely
	 * guarantees approximation quality. Otherwise, use a fast heuristic that has a higher chance of getting the
	 * estimate right the higher the number of samples
	 */
	ApproxBetweenness(const Graph& G, double epsilon=0.01, double delta=0.1, count diameterSamples=0);

	void run() override;

	/**
	 * @return number of samples taken in last run
	 */
	count numberOfSamples();


private:

	double epsilon;
	double delta;
	count r; // number of samples taken in last run
	count diameterSamples;
};

} /* namespace NetworKit */

#endif /* APPROXBETWEENNESS_H_ */
