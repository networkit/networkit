/*
 * GeneratorsBenchmark.cpp
 *
 *  Created on: May 29, 2013
 *      Author: forigem
 */

#include "GeneratorsBenchmark.h"

namespace NetworKit {

GeneratorsBenchmark::GeneratorsBenchmark() {
	// TODO Auto-generated constructor stub

}

GeneratorsBenchmark::~GeneratorsBenchmark() {
	// TODO Auto-generated destructor stub
}

TEST_F(GeneratorsBenchmark, benchmarkStaticBarabasiAlbertGenerator) {
	count k = 2;
	count nMax = 10e6;
	count n0 = 2;

	StaticBarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0);
	Graph G(0);
	EXPECT_TRUE(G.isEmpty());

	G = BarabasiAlbert.generate();
	EXPECT_FALSE(G.isEmpty());

	EXPECT_EQ(nMax, G.numberOfNodes());
	EXPECT_EQ( ((n0-1) + ((nMax - n0) * k)), G.numberOfEdges());



}

} /* namespace NetworKit */
