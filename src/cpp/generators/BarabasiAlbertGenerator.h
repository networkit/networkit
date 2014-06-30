/*
 * BarabasiAlbertGenerator.h
 *
 *  Created on: May 28, 2013
 *      Author: forigem
 */

#ifndef BarabasiAlbertGenerator_H_
#define BarabasiAlbertGenerator_H_

#include <set>

#include "StaticGraphGenerator.h"
#include "../auxiliary/ProgressMeter.h"

namespace NetworKit {

/**
 * @ingroup generators
 * Generates a scale-free graph using the Barabasi-Albert preferential attachment model.
 */
class BarabasiAlbertGenerator: public NetworKit::StaticGraphGenerator {
private:
	count k; //!< Attachments made per node
	count nMax; //!< The maximal number of nodes attached
	count n0; //!< The number of initial connected nodes
	count degreeSum; //!< Degree sum of the current graph


	Graph initializeGraph();

public:
	BarabasiAlbertGenerator();

	BarabasiAlbertGenerator(count k, count nMax, count n0 = 0);

	Graph generate() override;
};

} /* namespace NetworKit */
#endif /* BarabasiAlbertGenerator_H_ */
