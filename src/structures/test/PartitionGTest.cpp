/*
 * ClusteringAlgoGTest.cpp
 *
 *  Created on: 04.12.2013
 *      Author: Maximilian Vogel (uocvf@student.kit.edu)
 */

#include "PartitionGTest.h"

#include "../Partition.h"

#include <iostream>

#ifndef NOGTEST

namespace NetworKit {

/*TEST_F(PartitionGTest, test*) {
}*/

TEST_F(PartitionGTest, testAllToSingletons) {
	Partition p(10);
	p.allToSingletons();
	for (int i = 0; i < 10; i++) {
		EXPECT_TRUE(p[i] != p[i+1]);
	}
}

TEST_F(PartitionGTest, testAddToSubsetSuccess) {
	Partition p(10);
	p.addToSubset(8,2);
	p.addToSubset(3,2);
	EXPECT_TRUE(p[8] == p[3]);	
}

TEST_F(PartitionGTest, tryAddToSubsetWrongNodeID) {
	Partition p(10);
	p.addToSubset(30,15);
	// TODO: which behaviour is wanted? user/programmer has to take care of IDs?
}

TEST_F(PartitionGTest, tryAddToSubsetWrongSubsetID) {
	Partition p(10);
	p.addToSubset(15,8);
	EXPECT_TRUE(p[8] == none);
	// TODO: which behaviour is wanted? user/programmer has to take care of IDs with upperBound()?
}

TEST_F(PartitionGTest, testMoveToSubsetSuccess) {
	Partition p(10);
	p.allToSingletons();
	p.moveToSubset(p[3],8);
	EXPECT_TRUE(p[8] == p[3]);
}

TEST_F(PartitionGTest, tryMoveToSubsetWrongSubsetID) {
	Partition p(10);
	p.allToSingletons();
	p.moveToSubset(30,4);
	// TODO: which behaviour is wanted? user/programmer has to take care of IDs with upperBound()?
}

TEST_F(PartitionGTest, testMergeSubsets) {
	Partition p(10);
	p.allToSingletons();
	p.mergeSubsets(p[0],p[9]);
	p.mergeSubsets(p[1],p[8]);
	p.mergeSubsets(p[2],p[7]);
	p.mergeSubsets(p[0],p[1]);
	p.mergeSubsets(p[1],p[2]);
	EXPECT_TRUE(p[9] == p[7]);
}

TEST_F(PartitionGTest, testIsOnePartitionFalse) {
	Partition p(10);
	p.allToSingletons();
	EXPECT_FALSE(p.isOnePartition({}));
}

TEST_F(PartitionGTest, testIsOnePartitionFalse2) {
	Partition p(10);
	p.allToSingletons();
	p.forEntries([&](index e,index s){
		index n = (e==9)?0:e+1;
		p.mergeSubsets(p[e],p[n]);
	});
	p.toSingleton(9);
	EXPECT_FALSE(p.isOnePartition({}));
}

TEST_F(PartitionGTest, testIsOnePartitionTrue) {
	Partition p(10);
	p.allToSingletons();
	p.forEntries([&](index e,index s){
		p.moveToSubset(p[9],e);
	});
	EXPECT_TRUE(p.isOnePartition({}));
}

TEST_F(PartitionGTest, testIsSingletonPartitionTrue) {
	Partition p(10);
	p.allToSingletons();
	EXPECT_TRUE(p.isSingletonPartition({}));
}

TEST_F(PartitionGTest, testIsSingletonPartitionFalse) {
	Partition p(10);
	p.allToSingletons();
	p.mergeSubsets(p[0],p[9]);
	EXPECT_FALSE(p.isSingletonPartition({}));
}


TEST_F(PartitionGTest, testNumberOfSubsets) {
	count n = 6542;
	Partition p(n);
	p.allToSingletons();
	EXPECT_EQ(p.numberOfSubsets(),n);
}

TEST_F(PartitionGTest, testNumberOfSubsets2) {
	count n = 6542;
	Partition p(n);
	p.allToSingletons();
	p.parallelForEntries([&](index e, index s) {
		p.moveToSubset(p[0],e);
	});
	EXPECT_EQ(p.numberOfSubsets(),1);
}

TEST_F(PartitionGTest, testNumberOfSubsets3) {
	count n = 6542;
	Partition p(n);
	p.allToSingletons();
	for (int i = 0; i < 6542; i+=2) {
		p.mergeSubsets(p[i],p[i+1]);
	}
	EXPECT_EQ(p.numberOfSubsets(),3271);
}
TEST_F(PartitionGTest, testLowerBound) {
	count n = 6542;
	Partition p(n);
	p.allToSingletons();
	EXPECT_EQ(p.lowerBound(),0);
}

TEST_F(PartitionGTest, testUpperBound) {
	count n = 6542;
	Partition p(n);
	p.allToSingletons();
	// ID start with 1, plus 1 for upper bound
	count ub = n+2;
	EXPECT_EQ(p.upperBound(),ub);
}

TEST_F(PartitionGTest, testUpperBound2) {
	count n = 6542;
	Partition p(n);
	p.allToSingletons();
	for (int i = 0; i < 6542; i+=2) {
		p.mergeSubsets(p[i],p[i+1]);
	}
	// 6542 because of singletons
	// +(6542/2) merge operations
	// +1 for "hard" upper bound = 9815
	EXPECT_EQ(p.upperBound(),9815);
}

TEST_F(PartitionGTest, testContainsSuccessSingletons) {
	Partition p(10);
	p.allToSingletons();
	EXPECT_TRUE(p.contains(8));
}


TEST_F(PartitionGTest, testContainsSuccessAfterMerges) {
	Partition p(10);
	p.allToSingletons();
	p.forEntries([&](index e,index s){
		index n = (e==9)?0:e+1;
		p.mergeSubsets(p[e],p[n]);
	});
	p.toSingleton(9);
	EXPECT_TRUE(p.contains(8));
}

TEST_F(PartitionGTest, testContainsNodeOutOfRange) {
	Partition p(10);
	p.allToSingletons();
	EXPECT_FALSE(p.contains(15));
}

TEST_F(PartitionGTest, testContainsNodeWithoutSubset) {
	Partition p(10);
	for (int i = 0; i < 10; i++) {
		EXPECT_FALSE(p.contains(i));
	}
}

TEST_F(PartitionGTest, testInSameSubsetTrue) {
	Partition p(10);
	p.allToSingletons();
	p.mergeSubsets(p[0],p[9]);
	p.mergeSubsets(p[1],p[8]);
	p.mergeSubsets(p[2],p[7]);
	p.mergeSubsets(p[0],p[1]);
	p.mergeSubsets(p[1],p[2]);
	EXPECT_TRUE(p.inSameSubset(9,7));
}

TEST_F(PartitionGTest, testInSameSubsetFalse) {
	Partition p(10);
	p.allToSingletons();
	p.mergeSubsets(p[0],p[9]);
	p.mergeSubsets(p[1],p[8]);
	p.mergeSubsets(p[2],p[7]);
	p.mergeSubsets(p[0],p[1]);
	p.mergeSubsets(p[1],p[2]);
	EXPECT_FALSE(p.inSameSubset(9,4));
}

TEST_F(PartitionGTest, testCompact) {
	Partition p(10);
	p.allToSingletons();
	p.mergeSubsets(p[0],p[9]);
	p.mergeSubsets(p[1],p[8]);
	p.mergeSubsets(p[2],p[7]);
	p.mergeSubsets(p[0],p[1]);
	p.mergeSubsets(p[1],p[2]);
	p.compact();
	EXPECT_EQ(p.upperBound(),6);
}


} /* namespace NetworKit */

#endif /*NOGTEST */
