/*
 * CoarseningGTest.h
 *
 *  Created on: 20.12.2012
 *      Author: cls
 */

#ifndef COARSENINGGTEST_H_
#define COARSENINGGTEST_H_

#include <gtest/gtest.h>

#include "../../graph/Graph.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/ClusteringGenerator.h"
#include "../../coarsening/ClusterContracter.h"

namespace EnsembleClustering {

/**
 * googletest test fixture for the coarsening module.
 */
class CoarseningGTest: public testing::Test {

	// TODO: are constructor/destructor needed?

};

TEST_F(CoarseningGTest, testClusterContracter) {
	GraphGenerator graphGen;
	int64_t n = 10;
	Graph& G = graphGen.makeErdosRenyiGraph(n, 0.5);

	ClusteringGenerator clusteringGen;
	Clustering& singleton = clusteringGen.makeSingletonClustering(G);

	ClusterContracter contracter;
	Graph& Gcon = contracter.run(G, singleton);

	EXPECT_EQ(G.numberOfNodes(), Gcon.numberOfNodes())
			<< "graph contracted according to singleton clustering should have the same number of nodes as original";
	EXPECT_EQ(G.numberOfEdges(), Gcon.numberOfEdges())
			<< "graph contracted according to singletons clustering should have the same number of nodes as original";

	int k = 2; // number of clusters in random clustering
	Clustering& random = clusteringGen.makeRandomClustering(G, k);
	Graph& GConRand = contracter.run(G, random);

	EXPECT_EQ(k, GConRand.numberOfNodes())
			<< "graph contracted according to random clustering should have the same number of nodes as there are clusters.";

}

} /* namespace EnsembleClustering */
#endif /* COARSENINGGTEST_H_ */
