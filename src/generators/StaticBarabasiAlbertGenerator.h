/*
 * StaticBarabasiAlbertGenerator.h
 *
 *  Created on: May 28, 2013
 *      Author: forigem
 */

#ifndef STATICBARABASIALBERTGENERATOR_H_
#define STATICBARABASIALBERTGENERATOR_H_

#include <set>

#include "StaticGraphGenerator.h"
#include "../auxiliary/RandomInteger.h"
#include "../auxiliary/ProgressMeter.h"

namespace NetworKit {

class StaticBarabasiAlbertGenerator: public NetworKit::StaticGraphGenerator {
private:
	count k; //!< Attachments made per node
	count nMax; //!< The maximal number of nodes attached
	count n0; //!< The number of initial connected nodes
	count degreeSum; //!< Degree sum of the current graph
	Aux::RandomInteger randomInt;


	Graph initializeGraph();

public:
	StaticBarabasiAlbertGenerator(count k, count nMax, count n0 = -1);

	virtual ~StaticBarabasiAlbertGenerator();
	/**
	 * Generate a graph using the Barabasi-Albert preferential generator
	 */
	virtual Graph generate();
};

} /* namespace NetworKit */
#endif /* STATICBARABASIALBERTGENERATOR_H_ */
