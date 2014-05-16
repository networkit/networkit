/*
 * OverlapGTest.cpp
 *
 *  Created on: 21.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "OverlapGTest.h"
#include "../../auxiliary/Log.h"

namespace NetworKit {

TEST_F(OverlapGTest, testRegionGrowingOverlapperOnOneClustering) {
	GraphGenerator graphGen;
	int64_t n = 10;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusterGen;
	std::vector<Partition> clusterings;
	int z = 3; // number of clusterings
	for (int i = 0; i < z; ++i) {
 		clusterings.push_back(clusterGen.makeOneClustering(G));
	}
	DEBUG("end of loop");

	RegionGrowingOverlapper over;
	Partition core = over.run(G, clusterings);

	// test if core clustering is one-clustering
	node v = 1;
	index one = core.subsetOf(v);
	bool isOneClustering = true; //TODO replaces with function call?
	G.forNodes([&](node v) {
		index c = core.subsetOf(v);
		DEBUG("CLUSTER! c = " , c);
		isOneClustering = isOneClustering && (c == one);
	});

	EXPECT_TRUE(isOneClustering) << "overlap of multiple 1-clustering should be a 1-clustering";

}


TEST_F(OverlapGTest, testRegionGrowingOverlapperOnSingletonClustering) {
	GraphGenerator graphGen;
	int64_t n = 10;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusterGen;
	std::vector<Partition> clusterings;
	int z = 3; // number of clusterings
	for (int i = 0; i < z; ++i) {
		Partition zeta = clusterGen.makeSingletonClustering(G);
		clusterings.push_back(zeta);
	}

	RegionGrowingOverlapper over;
	Partition core = over.run(G, clusterings);

	// test if core clustering is singleton-clustering
	bool isSingleton = true; //TODO function call?
	G.forEdges([&](node u, node v) {
		isSingleton = isSingleton && (core.subsetOf(u) != core.subsetOf(v));
	});

	EXPECT_TRUE(isSingleton) << "overlap of multiple  singleton clusterings should be a singleton clustering";

}


TEST_F(OverlapGTest, testHashingOverlapperOnSingletonClusterings) {
	GraphGenerator graphGen;
	int64_t n = 10;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusterGen;
	std::vector<Partition> clusterings;
	int z = 3; // number of clusterings
	for (int i = 0; i < z; ++i) {
		clusterings.push_back(clusterGen.makeSingletonClustering(G));
	}

	HashingOverlapper over;
	Partition core = over.run(G, clusterings);

	// test if core clustering is singleton-clustering
	bool isSingleton = true;
	G.forEdges([&](node u, node v) {
		isSingleton = isSingleton && (core.subsetOf(u) != core.subsetOf(v));
	});

	EXPECT_TRUE(isSingleton) << "overlap of multiple  singleton clusterings should be a singleton clustering";
}


TEST_F(OverlapGTest, testHashingOverlapperOnOneClusterings) {
	GraphGenerator graphGen;
	int64_t n = 10;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusterGen;
	std::vector<Partition> clusterings;
	int z = 3; // number of clusterings
	for (int i = 0; i < z; ++i) {
 		clusterings.push_back(clusterGen.makeOneClustering(G));
	}
	DEBUG("end of loop");

	HashingOverlapper over;
	Partition core = over.run(G, clusterings);

	// test if core clustering is one-clustering
	node v = 1;
	index one = core.subsetOf(v);
	bool isOneClustering = true; //TODO function call?
	G.forNodes([&](node v) {
		index c = core.subsetOf(v);
		DEBUG("CLUSTER! c = " , c);
		isOneClustering = isOneClustering && (c == one);
	});

	EXPECT_TRUE(isOneClustering) << "overlap of multiple 1-clustering should be a 1-clustering";

}


TEST_F(OverlapGTest, testHashingOverlapperForCorrectness) {
	count n = 4;
	Graph G(n);

	Partition zeta(n);
	Partition eta(n);

	zeta.setUpperBound(2);
	zeta[0] = 0;
	zeta[1] = 0;
	zeta[2] = 1;
	zeta[3] = 1;

	eta.setUpperBound(2);
	eta[0] = 0;
	eta[1] = 1;
	eta[2] = 0;
	eta[3] = 1;

	std::vector<Partition> clusterings = {zeta, eta};
	HashingOverlapper overlapper;
	Partition overlap = overlapper.run(G, clusterings);

	INFO("overlap clustering number of clusters: ", overlap.numberOfSubsets());
	INFO("overlap clustering: ", overlap.getVector());
}




} /* namespace NetworKit */

#endif
