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
	bool batagelj; //!< Specifies whether to use batagelj's method or the original one

	Graph initializeGraph();

	/**
	 * Implementation of ALG 5 of https://kops.uni-konstanz.de/bitstream/handle/123456789/5799/random.pdf?sequence=1
	 * Running time is O(n+m)
	 * @return The generated graph
	 */
	Graph generateBatagelj();

public:
	BarabasiAlbertGenerator();

	/**
	 *
	 *
	 * @param k Number of attachments per node
	 * @param nMax Maximum number of nodes in the graph
	 * @param n0 Number of connected nodes to begin with
	 * @param batagelj Specifies whether to use batagelj's method or the original one; default: false
	 */
	BarabasiAlbertGenerator(count k, count nMax, count n0 = 0, bool batagelj=false);

	Graph generate() override;
};

} /* namespace NetworKit */
#endif /* BarabasiAlbertGenerator_H_ */
