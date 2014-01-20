/*
 * IndependentSetTest.cpp
 *
 *  Created on: 27.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST


#include "IndependentSetGTest.h"

namespace NetworKit {

IndependentSetGTest::IndependentSetGTest() {
	// TODO Auto-generated constructor stub

}

IndependentSetGTest::~IndependentSetGTest() {
	// TODO Auto-generated destructor stub
}


TEST_F(IndependentSetGTest, tryLuby) {
	count n = 500;
	GraphGenerator graphGen;
	Graph G = graphGen.makeRandomGraph(n, 0.001);
	INFO("G: " << G.toString());

	Luby luby;
	std::vector<bool> I = luby.run(G);

	EXPECT_TRUE(luby.isIndependentSet(I, G)) << "result must be an independent set";

	count size = 0;
	for (bool x : I) {
		if (x) {
			size += 1;
		}
	}
	INFO("independent set size: " << size << "/" << n);
}

TEST_F(IndependentSetGTest, tryLubyWithSelfLoops) {
	count n = 500;
	GraphGenerator graphGen;
	Graph G = graphGen.makeRandomGraph(n, 0.001);
	G.forNodes([&](node u){
		G.addEdge(u,u);
	});

	INFO("G: " << G.toString());

	Luby luby;
	std::vector<bool> I = luby.run(G);

	EXPECT_TRUE(luby.isIndependentSet(I, G)) << "result must be an independent set";

	count size = 0;
	for (bool x : I) {
		if (x) {
			size += 1;
		}
	}
	INFO("independent set size: " << size << "/" << n);
}

} /* namespace NetworKit */

#endif /*NOGTEST */

