/*
 * GeneratorsBenchmark.cpp
 *
 *  Created on: May 29, 2013
 *      Author: forigem
 */

#ifndef NOGTEST

#include "GeneratorsBenchmark.h"
#include "../../auxiliary/Log.h"

#include "../HyperbolicGenerator.h"
#include "../DynamicHyperbolicGenerator.h"

namespace NetworKit {

GeneratorsBenchmark::GeneratorsBenchmark() {
	// TODO Auto-generated constructor stub

}

GeneratorsBenchmark::~GeneratorsBenchmark() {
	// TODO Auto-generated destructor stub
}

TEST_F(GeneratorsBenchmark, benchmarkBarabasiAlbertGenerator) {
	count k = 2;
	count nMax = 100000;
	count n0 = 2;

	BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0);
	Graph G(0);
	EXPECT_TRUE(G.isEmpty());

	G = BarabasiAlbert.generate();
	EXPECT_FALSE(G.isEmpty());

	EXPECT_EQ(nMax, G.numberOfNodes());
	EXPECT_EQ( ((n0-1) + ((nMax - n0) * k)), G.numberOfEdges());

}

TEST_F(GeneratorsBenchmark, benchmarkHyperbolicGenerator) {
	count n = 500000;
	HyperbolicGenerator gen(n,1,1);
	Graph G = gen.generate();
	EXPECT_EQ(G.numberOfNodes(), n);
}

TEST_F(GeneratorsBenchmark, benchmarkDynamicHyperbolicGenerator) {
	count n = 10000;
	count nSteps = 100;
	//(count n, double initialFactor = 1, double alpha = 1, double stretch = 1, double moveEachStep = 0, double factorgrowth = 0, double moveDistance = 0);
	DynamicHyperbolicGenerator dyngen(n,0,1,1,0,0.01,0);
	dyngen.generate(nSteps);

}

} /* namespace NetworKit */

#endif /*NOGTEST */
