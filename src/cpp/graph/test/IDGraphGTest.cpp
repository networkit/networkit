/*
 * IDGraphGTest.cpp
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#include "IDGraphGTest.h"

#include "../Graph.h"
#include "../DirectedGraph.h"
#include "../BasicGraph.h"

namespace NetworKit {

using testing::Types;
typedef Types<DirectedGraph, DirectedGraph_T> dgraphImplementations;
TYPED_TEST_CASE(IDGraphGTest, dgraphImplementations);

template <typename T>
void IDGraphGTest<T>::SetUp() {
	/*
	 *    0
	 *   . \
	 *  /   \
	 * /     .
	 * 1 <-- 2
	 * ^ \  .|
	 * |  \/ |
	 * | / \ |
	 * |/   ..
	 * 3 <-- 4
	 *
	 * move you pen from node to node:
	 * 3 -> 1 -> 0 -> 2 -> 1 -> 4 -> 3 -> 2 -> 4
	 */
	Ghouse = T(5);
	houseEdgesOut = {
		{0, 2},
		{1, 0},
		{1, 4},
		{2, 1},
		{2, 4},
		{3, 1},
		{3, 2},
		{4, 3}
	};
	for (auto& e : houseEdgesOut) {
		Ghouse.addEdge(e.first, e.second);
	}
	m_house = Ghouse.numberOfEdges();
}

TYPED_TEST(IDGraphGTest, degreeInOut) {
	std::vector<count> inDegrees(this->Ghouse.upperNodeIdBound(), 0);
	std::vector<count> outDegrees(this->Ghouse.upperNodeIdBound(), 0);
	for (auto& e : this->houseEdgesOut) {
		outDegrees[e.first]++;
		inDegrees[e.second]++;
	}

	this->Ghouse.forNodes([&](node v) {
		EXPECT_EQ(inDegrees[v], this->Ghouse.degreeIn(v));
		EXPECT_EQ(outDegrees[v], this->Ghouse.degreeOut(v));
		EXPECT_EQ(outDegrees[v], this->Ghouse.degree(v));
	});
}

// TYPED_TEST(IDGraphGTest, forOutEdgesOf) {
// 	count m = 0;
// 	std::vector<bool> visited(this->m_house, false);

// 	this->Ghouse.forNodes([&](node u) {
// 		this->Ghouse.forOutEdgesOf(u, [&](node v, node w) {
// 			// edges should be v to w, so if we iterate over edges from u, u should be equal v
// 			EXPECT_EQ(u, v);
			
// 			auto e = std::make_pair(v, w);
// 			// find edge
// 			auto it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);

// 			EXPECT_FALSE(it == this->houseEdgesOut.end()); // check if edge is allowed to exists
		
// 			// find index in edge array
// 			int i = std::distance(this->houseEdgesOut.begin(), it);
// 			EXPECT_FALSE(visited[i]); // make sure edge was not visited before (would be visited twice)
			
// 			// mark edge as visited
// 			visited[i] = true;
// 			m++;
// 		});
// 	});

// 	EXPECT_EQ(this->m_house, m);
// 	for (auto b : visited) {
// 		EXPECT_TRUE(b);
// 	}
// }

// TYPED_TEST(IDGraphGTest, forInEdgesOf) {
// 	// very similar to forOutEdgesOf ...

// 	count m = 0;
// 	std::vector<bool> visited(this->m_house, false);

// 	// NEXT 3 LINES ARE DIFFERENT
// 	this->Ghouse.forNodes([&](node u) {
// 		this->Ghouse.forInEdgesOf(u, [&](node v, node w) {
// 			// edges should be v to w, so if we iterate over edges from u, u should be equal w
// 			EXPECT_EQ(u, w);
			
// 			auto e = std::make_pair(v, w);
// 			// find edge
// 			auto it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
// 			EXPECT_FALSE(it == this->houseEdgesOut.end()); // check if edge is allowed to exists
		
// 			// find index in edge array
// 			int i = std::distance(this->houseEdgesOut.begin(), it);
// 			EXPECT_FALSE(visited[i]); // make sure edge was not visited before (would be visited twice)
			
// 			// mark edge as visited
// 			visited[i] = true;
// 			m++;
// 		});
// 	});

// 	EXPECT_EQ(this->m_house, m);
// 	for (auto b : visited) {
// 		EXPECT_TRUE(b);
// 	}
// }

TYPED_TEST(IDGraphGTest, BFSfrom) {
	std::vector<count> visitedOrder(5, none);
	index i = 0;
	this->Ghouse.BFSfrom(3, [&](node v) {
		EXPECT_EQ(none, visitedOrder[v]); // visit every node once
		visitedOrder[v] = i++;
	});
	// have we visited all nodes
	for (count l : visitedOrder) {
		EXPECT_TRUE(l != none);
	}

	// root on level 0
	EXPECT_EQ(0u, visitedOrder[3]);

	// level 1
	EXPECT_TRUE( (visitedOrder[1] == 1) ^ (visitedOrder[1] == 2) );
	EXPECT_TRUE( (visitedOrder[2] == 1) ^ (visitedOrder[2] == 2) );

	// level 2
	EXPECT_TRUE( (visitedOrder[0] == 3) ^ (visitedOrder[0] == 4) );
	EXPECT_TRUE( (visitedOrder[4] == 3) ^ (visitedOrder[4] == 4) );
}

TYPED_TEST(IDGraphGTest, DFSfrom) {
	std::vector<count> visitedOrder(5, none);
	index i = 0;
	this->Ghouse.DFSfrom(3, [&](node v) {
		EXPECT_EQ(none, visitedOrder[v]); // visit every node once
		visitedOrder[v] = i++;
	});

	// have we visited all nodes
	for (count l : visitedOrder) {
		EXPECT_TRUE(l != none);
	}

	// root on level 0
	EXPECT_EQ(0u, visitedOrder[3]);

	// level 1
	EXPECT_TRUE( (visitedOrder[1] == 1) ^ (visitedOrder[2] == 1) );

	// level 2
	EXPECT_TRUE( (visitedOrder[0] == 2) ^ (visitedOrder[1] == 2) ^ (visitedOrder[4] == 2) );

	// level 3
	EXPECT_TRUE( (visitedOrder[2] == 3) ^ (visitedOrder[0] == 3) ^ (visitedOrder[4] == 3) ^ (visitedOrder[1] == 3) );

	// level 4
	EXPECT_TRUE( (visitedOrder[2] == 4) ^ (visitedOrder[4] == 4) ^ (visitedOrder[0] == 4) );
}

} /* namespace NetworKit */

#endif /*NOGTEST */
