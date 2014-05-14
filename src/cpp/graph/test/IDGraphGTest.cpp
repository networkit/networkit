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

namespace NetworKit {

using testing::Types;
typedef Types<DirectedGraph> dgraphImplementations;
TYPED_TEST_CASE(IDGraphGTest, dgraphImplementations);

template <typename T>
void IDGraphGTest<T>::SetUp() {
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
		ASSERT_EQ(inDegrees[v], this->Ghouse.degreeIn(v));
		ASSERT_EQ(outDegrees[v], this->Ghouse.degreeOut(v));
		// ASSERT_EQ(outDegrees[v], this->Ghouse.degree(v));
	});
}

TYPED_TEST(IDGraphGTest, forOutEdgesOf) {
	count m = 0;
	std::vector<bool> visited(this->m_house, false);

	this->Ghouse.forNodes([&](node u) {
		this->Ghouse.forOutEdgesOf(u, [&](node v, node w) {
			// edges should be v to w, so if we iterate over edges from u, u should be equal v
			ASSERT_EQ(u, v);
			
			auto e = std::make_pair(v, w);
			// find edge
			auto it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
			ASSERT_FALSE(it == this->houseEdgesOut.end()); // check if edge is allowed to exists
		
			// find index in edge array
			int i = std::distance(this->houseEdgesOut.begin(), it);
			ASSERT_FALSE(visited[i]); // make sure edge was not visited before (would be visited twice)
			
			// mark edge as visited
			visited[i] = true;
			m++;
		});
	});

	ASSERT_EQ(this->m_house, m);
	for (auto b : visited) {
		ASSERT_TRUE(b);
	}
}

TYPED_TEST(IDGraphGTest, forInEdgesOf) {
	// very similar to forOutEdgesOf ...

	count m = 0;
	std::vector<bool> visited(this->m_house, false);

	// NEXT 3 LINES ARE DIFFERENT
	this->Ghouse.forNodes([&](node u) {
		this->Ghouse.forInEdgesOf(u, [&](node v, node w) {
			// edges should be v to w, so if we iterate over edges from u, u should be equal w
			ASSERT_EQ(u, w);
			
			auto e = std::make_pair(v, w);
			// find edge
			auto it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
			ASSERT_FALSE(it == this->houseEdgesOut.end()); // check if edge is allowed to exists
		
			// find index in edge array
			int i = std::distance(this->houseEdgesOut.begin(), it);
			ASSERT_FALSE(visited[i]); // make sure edge was not visited before (would be visited twice)
			
			// mark edge as visited
			visited[i] = true;
			m++;
		});
	});

	ASSERT_EQ(this->m_house, m);
	for (auto b : visited) {
		ASSERT_TRUE(b);
	}
}

} /* namespace NetworKit */

#endif /*NOGTEST */
