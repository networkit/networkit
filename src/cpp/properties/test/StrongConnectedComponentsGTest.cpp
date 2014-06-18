/*
 * StrongConnectedComponentsGTest.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */
#ifndef NOGTEST

#include "StrongConnectedComponentsGTest.h"
#include "../StrongConnectedComponents.h"

namespace NetworKit {

StrongConnectedComponentsGTest::StrongConnectedComponentsGTest() {

}

StrongConnectedComponentsGTest::~StrongConnectedComponentsGTest() {

}

 TEST_F(StrongConnectedComponentsGTest, testWikiExample) {
 	count n = 8;
 	count m = 14;
 	Graph G(n, false, true);

 	G.addEdge(0, 4);
 	G.addEdge(1, 0);
 	G.addEdge(2, 1);
 	G.addEdge(2, 3);
 	G.addEdge(3, 2);
 	G.addEdge(4, 1);
 	G.addEdge(5, 1);
 	G.addEdge(5, 4);
 	G.addEdge(5, 6);
 	G.addEdge(6, 2);
 	G.addEdge(6, 5);
 	G.addEdge(7, 3);
 	G.addEdge(7, 6);
 	G.addEdge(7, 7);

 	ASSERT_EQ(n, G.numberOfNodes());
 	ASSERT_EQ(m, G.numberOfEdges());

 	count z = G.upperNodeIdBound();
	Partition p_expected(z);
	p_expected.allToSingletons();
	p_expected[0] = 0;
	p_expected[1] = 0;
	p_expected[2] = 1;
	p_expected[3] = 1;
	p_expected[4] = 0;
	p_expected[5] = 2;
	p_expected[6] = 2;
	p_expected[7] = 3;
	p_expected.compact();

	StrongConnectedComponents scc(G);
	scc.run();
	Partition p_actual = scc.getPartition();
	p_actual.compact();

	comparePartitions(p_expected, p_actual);
}

void StrongConnectedComponentsGTest::comparePartitions(const Partition& p1, const Partition& p2) const {
	std::vector<index> partitionIdMap(p1.upperBound(), none);
	ASSERT_EQ(p1.numberOfElements(), p2.numberOfElements());
	ASSERT_EQ(p1.numberOfSubsets(), p2.numberOfSubsets());

	p1.forEntries([&](node v, index p) {
		if (partitionIdMap[p] == none) {
			partitionIdMap[p] = p2.subsetOf(v);
		}
		index p_mapped = partitionIdMap[p];
		ASSERT_EQ(p_mapped, p);
	});
}

} /* namespace NetworKit */

#endif /*NOGTEST */

