/*
 * CoarseningGTest.cpp
 *
 *  Created on: 20.12.2012
 *      Author: cls
 */

#include "CoarseningGTest.h"

namespace EnsembleClustering {

TEST_F(CoarseningGTest, testClusterContracter) {
	GraphGenerator graphGen;
	int64_t n = 100;
	Graph G = graphGen.makeErdosRenyiGraph(n, 0.5);

	ClusteringGenerator clusteringGen;
	Clustering singleton = clusteringGen.makeSingletonClustering(G);


	ClusterContracter contracter;
	GraphContraction conSingleton = contracter.run(G, singleton);
	Graph Gcon = conSingleton.getCoarseGraph();

	EXPECT_EQ(G.numberOfNodes(), Gcon.numberOfNodes())
			<< "graph contracted according to singleton clustering should have the same number of nodes as original";
	EXPECT_EQ(G.numberOfEdges(), Gcon.numberOfEdges())
			<< "graph contracted according to singletons clustering should have the same number of nodes as original";

	int k = 2; // number of clusters in random clustering
	Clustering random = clusteringGen.makeRandomClustering(G, k);
	GraphContraction conRand = contracter.run(G, random);
	Graph GconRand = conRand.getCoarseGraph();

	EXPECT_EQ(k, GconRand.numberOfNodes())
			<< "graph contracted according to random clustering should have the same number of nodes as there are clusters.";

}


TEST_F(CoarseningGTest, testClusteringProjector) {
	GraphGenerator graphGen;
	int64_t n = 100;
	Graph G = graphGen.makeErdosRenyiGraph(n, 0.5);

	// make random clustering with k clusters
	int k = 10;
	ClusteringGenerator clusteringGen;
	Clustering zetaRand = clusteringGen.makeRandomClustering(G, k);

	ClusterContracter contracter;
	GraphContraction contraction = contracter.run(G, zetaRand);
	EnsembleClustering::Clustering zetaOne = clusteringGen.makeOneClustering(contraction.getCoarseGraph());

}



} /* namespace EnsembleClustering */
