/*
 * OverlapGTest.h
 *
 *  Created on: 21.12.2012
 *      Author: cls
 */

#ifndef OVERLAPGTEST_H_
#define OVERLAPGTEST_H_

#include <gtest/gtest.h>

#include "../RegionGrowingOverlapper.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/ClusteringGenerator.h"


namespace EnsembleClustering {

class OverlapGTest: public testing::Test {

};


TEST_F(OverlapGTest, testRegionGrowingOverlapperOnOneClustering) {
	GraphGenerator graphGen;
	int64_t n = 10;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusterGen;
	std::vector<Clustering> clusterings;
	int z = 3; // number of clusterings
	for (int i = 0; i < z; ++i) {
 		clusterings.push_back(clusterGen.makeOneClustering(G));
	}
	DEBUG("end of loop");

	RegionGrowingOverlapper over;
	Clustering core = over.run(G, clusterings);

	// test if core clustering is one-clustering
	node v = 1;
	cluster one = core.clusterOf(v);
	bool isOneClustering = true;
	G.forallNodes([&](node v) {
		cluster c = core.clusterOf(v);
		DEBUG("CLUSTER! c = " << c);
		isOneClustering = isOneClustering && (c == one);
	});

	EXPECT_TRUE(isOneClustering) << "overlap of multiple 1-clustering should be a 1-clustering";

}


TEST_F(OverlapGTest, testRegionGrowingOverlapperOnSingletonClustering) {
	GraphGenerator graphGen;
	int64_t n = 10;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusterGen;
	std::vector<Clustering> clusterings;
	int z = 3; // number of clusterings
	for (int i = 0; i < z; ++i) {
		Clustering zeta = clusterGen.makeSingletonClustering(G);
		clusterings.push_back(zeta);
	}

	RegionGrowingOverlapper over;
	Clustering core = over.run(G, clusterings);

	// test if core clustering is singleton-clustering
	bool isSingleton = true;
	G.forallEdges([&](node u, node v) {
		isSingleton = isSingleton && (core.clusterOf(u) != core.clusterOf(v));
	});

	EXPECT_TRUE(isSingleton) << "overlap of multiple  singleton clusterings should be a singleton clustering";

}

} /* namespace EnsembleClustering */
#endif /* OVERLAPGTEST_H_ */
