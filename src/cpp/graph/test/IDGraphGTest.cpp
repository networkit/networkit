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

} /* namespace NetworKit */

#endif /*NOGTEST */
