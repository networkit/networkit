/*
 * ClusteringTest.cpp
 *
 *  Created on: 12.12.2012
 *      Author: cls
 */

#include "ClusteringGTest.h"

namespace EnsembleClustering {



TEST_F(ClusteringGTest, testModularity) {
	GraphGenerator graphGenerator;

	int n = 100;

	DEBUG("testing modularity on clustering of complete graph with " << n << " nodes");


	Graph G = graphGenerator.makeCompleteGraph(n);
	DEBUG("total edge weight: " << G.totalEdgeWeight());

	ClusteringGenerator clusteringGenerator;

	Clustering singleton = clusteringGenerator.makeSingletonClustering(G);
	Clustering one = clusteringGenerator.makeOneClustering(G);

	Modularity modularity;

	DEBUG("calculating modularity for singleton clustering");
	double modSingleton = modularity.getQuality(singleton, G);

	DEBUG("calculating modularity for 1-clustering");
	double modOne = modularity.getQuality(one, G);

	DEBUG("mod(singleton-clustering) = " << modSingleton);
	DEBUG("mod(1-clustering) = " << modOne);


	EXPECT_EQ(0.0, modOne) << "1-clustering should have modularity of 0.0";
	EXPECT_GE(0.0, modSingleton) << "singleton clustering should have modularity less than 0.0";

}

TEST_F(ClusteringGTest, testCoverage) {
	GraphGenerator graphGenerator;

	int n = 100;

	DEBUG("testing coverage on clustering of complete graph with " << n << " nodes");


	Graph G = graphGenerator.makeCompleteGraph(n);

	ClusteringGenerator clusteringGenerator;

	Clustering singleton = clusteringGenerator.makeSingletonClustering(G);
	Clustering one = clusteringGenerator.makeOneClustering(G);

	Coverage coverage;

	DEBUG("calculating coverage for singleton clustering");
	double covSingleton = coverage.getQuality(singleton, G);

	DEBUG("calculating coverage for 1-clustering");
	double covOne = coverage.getQuality(one, G);

	DEBUG("mod(singleton-clustering) = " << covSingleton);
	DEBUG("mod(1-clustering) = " << covOne);


	EXPECT_EQ(1.0, covOne) << "1-clustering should have coverage of 1.0";
	EXPECT_GE(0.0, covSingleton) << "singleton clustering should have coverage of 0.0";

}


TEST_F(ClusteringGTest, testClusteringEquality) {
	int n = 100;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Clustering one1 = clusteringGen.makeOneClustering(G);
	Clustering one2 = clusteringGen.makeOneClustering(G);

	EXPECT_TRUE(one1.equals(one2, G)) << "two 1-clusterings of G should be equal";

	Clustering singleton1 = clusteringGen.makeSingletonClustering(G);
	Clustering singleton2 = clusteringGen.makeSingletonClustering(G);

	EXPECT_TRUE(singleton1.equals(singleton2, G)) << "two singleton clusterings of G should be equal";

}



TEST_F(ClusteringGTest, testJaccardMeasure) {
	int n = 100;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Clustering singleton = clusteringGen.makeSingletonClustering(G);
	Clustering random = clusteringGen.makeRandomClustering(G, 10);

	JaccardMeasure jaccard;
	double j = jaccard.getDissimilarity(G, singleton, random);

	EXPECT_EQ(1.0, j) << "The singleton clustering and any other clustering compare with a dissimilarity of 1.0";

}


TEST_F(ClusteringGTest, testRandMeasure) {
	int n = 100;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Clustering one1 = clusteringGen.makeOneClustering(G);
	Clustering one2 = clusteringGen.makeOneClustering(G);

	RandMeasure rand;
	double r = rand.getDissimilarity(G, one1, one2);

	EXPECT_EQ(0.0, r) << "Identical clusterings should compare with a dissimilarity of 0.0";

}


TEST_F(ClusteringGTest, testCompact) {
	int n = 50;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Clustering clustering = clusteringGen.makeRandomClustering(G, 2*n);
	clustering.compact();

	EXPECT_LE(clustering.upperBound(), n);
}



} /* namespace EnsembleClustering */
