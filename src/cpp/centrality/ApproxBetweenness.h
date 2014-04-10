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

	ApproxBetweenness(const Graph& G, bool normalized=false, double epsilon=0, double delta=0);

	void run() override;


private:

	double epsilon;
	double delta;
};

} /* namespace NetworKit */

#endif /* APPROXBETWEENNESS_H_ */
