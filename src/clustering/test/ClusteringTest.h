/*
 * ClusteringTest.h
 *
 *  Created on: 12.12.2012
 *      Author: cls
 */

#ifndef CLUSTERINGTEST_H_
#define CLUSTERINGTEST_H_

#include <gtest/gtest.h>


#include "../../aux/log.h"
#include "../Clustering.h"
#include "../Modularity.h"
#include "../ClusteringGenerator.h"
#include "../../graph/GraphGenerator.h"
#include "../LabelPropagation.h"

namespace EnsembleClustering {

class ClusteringTest: public testing::Test {


};



TEST_F(ClusteringTest, testModularity) {
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


TEST_F(ClusteringTest, testLabelPropagationOnUniformGraph) {
	GraphGenerator graphGenerator;
	int n = 100;
	Graph G = graphGenerator.makeErdosRenyiGraph(n, 0.2);

	LabelPropagation lp;
	Clustering zeta = lp.run(G);

	EXPECT_TRUE(zeta.isProper(G)) << "the resulting partition should be a proper clustering";

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	DEBUG("modularity produced by LabelPropagation: " << mod);
	EXPECT_GE(1.0, mod) << "valid modularity values are in [-0.5, 1]";
	EXPECT_LE(-0.5, mod) << "valid modularity values are in [-0.5, 1]";
}


TEST_F(ClusteringTest, testLabelPropagationOnClusteredGraph) {
	GraphGenerator graphGenerator;
	int n = 100;
	int k = 3; // number of clusters
	Graph G = graphGenerator.makeClusteredRandomGraph(n, k, 1.0, 0.001);

	LabelPropagation lp;
	Clustering zeta = lp.run(G);

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	DEBUG("modularity produced by LabelPropagation: " << mod);

	EXPECT_TRUE(zeta.isProper(G)) << "the resulting partition should be a proper clustering";
	EXPECT_EQ(k, zeta.numberOfClusters()) << " " << k << " clusters are easy to detect";

}


} /* namespace EnsembleClustering */
#endif /* CLUSTERINGTEST_H_ */
