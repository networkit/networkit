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




} /* namespace EnsembleClustering */
