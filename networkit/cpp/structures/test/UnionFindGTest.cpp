/*
 * PartitionGTest.cpp
 *
 *  Created on: 04.12.2013
 *      Author: Maximilian Vogel (uocvf@student.kit.edu)
 */

#include <iostream>

#include "UnionFindGTest.h"

#include "../UnionFind.h"

#ifndef NOGTEST

namespace NetworKit {

TEST_F(UnionFindGTest, testAllToSingletons) {
	UnionFind p(10);
	p.allToSingletons();
	for (int i = 0; i < 10; i++) {
		EXPECT_TRUE(p.find(i) != p.find((i+1)%10));
	}
}

TEST_F(UnionFindGTest, testMergeSimple) {
	UnionFind p(10);
	p.merge(3,8);
	EXPECT_EQ(p.find(8),p.find(3));	
}


TEST_F(UnionFindGTest, testMergeSubsets) {
	UnionFind p(10);
	p.allToSingletons();
	p.merge(0,9);
	p.merge(1,8);
	p.merge(2,7);
	p.merge(0,1);
	p.merge(1,2);
	EXPECT_EQ(p.find(9),p.find(7));
}

} /* namespace NetworKit */

#endif /*NOGTEST */
