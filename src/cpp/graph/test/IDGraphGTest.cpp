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

template <typename T>
IDGraphGTest<T>::IDGraphGTest() {
	// TODO Auto-generated constructor stub

}

template <typename T>
IDGraphGTest<T>::~IDGraphGTest() {
	// TODO Auto-generated destructor stub
}

using testing::Types;
typedef Types<DirectedGraph> dgraphImplementations;

TYPED_TEST_CASE(IDGraphGTest, dgraphImplementations);

TYPED_TEST(IDGraphGTest, degreeInOut) {
	TypeParam G(5);

	G.addEdge(0, 1);
	G.addEdge(1, 3);
	G.addEdge(3, 1);
	G.addEdge(4, 3);

	ASSERT_EQ(0, G.degreeIn(0));
	ASSERT_EQ(2, G.degreeIn(1));
	ASSERT_EQ(0, G.degreeIn(2));
	ASSERT_EQ(2, G.degreeIn(3));
	ASSERT_EQ(0, G.degreeIn(4));
	ASSERT_EQ(1, G.degreeOut(0));
	ASSERT_EQ(1, G.degreeOut(1));
	ASSERT_EQ(0, G.degreeOut(2));
	ASSERT_EQ(1, G.degreeOut(3));
	ASSERT_EQ(1, G.degreeOut(4));

	G.addEdge(3, 3);
	G.removeEdge(0, 1);

	ASSERT_EQ(0, G.degreeIn(0));
	ASSERT_EQ(1, G.degreeIn(1));
	ASSERT_EQ(0, G.degreeIn(2));
	ASSERT_EQ(3, G.degreeIn(3));
	ASSERT_EQ(0, G.degreeIn(4));
	ASSERT_EQ(0, G.degreeOut(0));
	ASSERT_EQ(1, G.degreeOut(1));
	ASSERT_EQ(0, G.degreeOut(2));
	ASSERT_EQ(2, G.degreeOut(3));
	ASSERT_EQ(1, G.degreeOut(4));
}


TYPED_TEST(IDGraphGTest, IterateNeighborMethodsDirected) {

	TypeParam G(5);
	std::vector<std::vector<node>> visited(5, std::vector<node>(5));
	std::vector<node> visitedNode(5);
	
	
	G.addEdge(0, 2);
	G.addEdge(1, 0);
	G.addEdge(1, 4);	
	G.addEdge(2, 1);
	G.addEdge(2, 4);
	G.addEdge(3, 1);
	G.addEdge(3, 2);
	G.addEdge(4, 3);

	
	G.forOutEdgesOf(0,[&](node u, node v){
		
		visited[u][v]=1;
	});


	ASSERT_EQ(1,visited[0][2]);

	G.forOutEdgesOf(0,[&](node u, node v){
		
		visited[u][v]=0;
	});


	G.forInEdgesOf(4, [&](node v, node u){

		visited[u][v] = 1;
	});

	ASSERT_EQ(1,visited[2][4]);
	ASSERT_EQ(1,visited[1][4]);

	G.forOutNeighborsOf(3, [&](node v){
		
		visitedNode[v]= 1;
		
	});
	
	G.forInNeighborsOf(2, [&](node v){
		
		visitedNode[v]= 1;
		
	});

	ASSERT_EQ(1, visitedNode[1]);
	ASSERT_EQ(1, visitedNode[2]);
	ASSERT_EQ(1, visitedNode[3]);
	ASSERT_EQ(1, visitedNode[0]);
}

} /* namespace NetworKit */

#endif /*NOGTEST */
