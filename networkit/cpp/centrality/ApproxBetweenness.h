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
	 * The values are normalized by default. The run() method takes O(m) time per sample, where  m is
	 * the number of edges of the graph. The number of samples is proportional to universalConstant/epsilon^2.
	 *
	 * @param	G			the graph
	 * @param	epsilon		maximum additive error
	 * @param	delta		probability that the values are within the error guarantee
	 * @param	diameterSamples		if 0, use the possibly slow estimation of the vertex diameter which definitely
	 * guarantees approximation quality. Otherwise, use a fast heuristic that has a higher chance of getting the
	 * estimate right the higher the number of samples (note: there is no
	 * approximation guarantee if using the heuristic)
	 * @param   universalConstant   the universal constant to be used in
	 * computing the sample size. It is 1 by default. Some references suggest
	 * using 0.5, but there is no guarantee in this case.
	 */
	ApproxBetweenness(const Graph& G, const double epsilon=0.01, const double delta=0.1, const count diameterSamples=0, const double universalConstant=1.0);

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
	double universalConstant;
};

} /* namespace NetworKit */

#endif /* APPROXBETWEENNESS_H_ */
