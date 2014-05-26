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
 * Approximation of betweenness centrality according to algorithm described in
 * 	Matteo Riondato and Evgenios M. Kornaropoulos: Fast Approximation of Betweenness Centrality through Sampling
 */
class ApproxBetweenness: public NetworKit::Centrality {

public:

	/** 
	 * The algorithm approximates the betweenness of all vertices so that the scores are
	 * within an additive error epsilon with probability at least (1- delta).
	 * The values are normalized by default.
	 * 
	 * @param	epsilon		maximum additive error
	 * @param	delta		probability that the values are within the error guarantee
	 */
	ApproxBetweenness(const Graph& G, double epsilon=0.01, double delta=0.1);

	void run() override;


private:

	double epsilon;
	double delta;
};

} /* namespace NetworKit */

#endif /* APPROXBETWEENNESS_H_ */
