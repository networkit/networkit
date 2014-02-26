/*
 * ClusteringAlgoGTest.cpp
 *
 *  Created on: 10.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ClusteringAlgoGTest.h"

#include "../PLP.h"
#include "../PLM.h"
#include "../CNM.h"
#include "../ParallelAgglomerativeClusterer.h"
#include "../../clustering/Modularity.h"
#include "../../graph/GraphGenerator.h"
#include "../../clustering/ClusteringGenerator.h"
#include "../../io/METISGraphReader.h"
#include "../EPP.h"
#include "../../overlap/HashingOverlapper.h"
#include "../EPPFactory.h"
#include "../CommunityGraph.h"
#include "../PLM2.h"
#include "../../clustering/GraphClusteringTools.h"



#ifndef NOGTEST

namespace NetworKit {

TEST_F(ClusteringAlgoGTest, testEnsemblePreprocessing) {
	count n = 1000;
	count k = 10;
	double pin = 1.0;
	double pout = 0.0;

	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	EPP ensemble;

	count b = 4;
	for (count i = 0; i < b; ++i) {
		ensemble.addBaseClusterer(*(new PLP()));
	}
	ensemble.setFinalClusterer(*(new PLM()));
	ensemble.setOverlapper(*(new HashingOverlapper));

	Partition zeta = ensemble.run(G);

	INFO("number of clusters:" , zeta.numberOfSubsets());

	Modularity modularity;
	INFO("modularity: " , modularity.getQuality(zeta, G));



}



TEST_F(ClusteringAlgoGTest, testLabelPropagationOnUniformGraph) {
	GraphGenerator graphGenerator;
	int n = 100;
	Graph G = graphGenerator.makeErdosRenyiGraph(n, 0.2);

	PLP lp;
	Partition zeta = lp.run(G);

	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta)) << "the resulting partition should be a proper clustering";

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	DEBUG("modularity produced by LabelPropagation: " , mod);
	EXPECT_GE(1.0, mod) << "valid modularity values are in [-0.5, 1]";
	EXPECT_LE(-0.5, mod) << "valid modularity values are in [-0.5, 1]";
}


TEST_F(ClusteringAlgoGTest, testLabelPropagationOnClusteredGraph_ForNumberOfClusters) {
	GraphGenerator graphGenerator;
	int64_t n = 100;
	count k = 3; // number of clusters
	Graph G = graphGenerator.makeClusteredRandomGraph(n, k, 1.0, 0.001);

	PLP lp;
	Partition zeta = lp.run(G);

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	DEBUG("modularity produced by LabelPropagation: " , mod);

	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta)) << "the resulting partition should be a proper clustering";
	EXPECT_EQ(k, zeta.numberOfSubsets()) << " " << k << " clusters are easy to detect";

}


TEST_F(ClusteringAlgoGTest, testLabelPropagationOnClusteredGraph_ForEquality) {
	int64_t n = 100;

	GraphGenerator graphGen;
	Graph Gtrash = graphGen.makeCompleteGraph(n);

	count k = 3; // number of clusters
	ClusteringGenerator clusteringGen;
	Partition reference = clusteringGen.makeRandomClustering(Gtrash, k);
	assert (reference.numberOfSubsets() == k);

	Graph G = graphGen.makeClusteredRandomGraph(reference, 1.0, 0.0);	// LabelPropagation is very bad at discerning clusters and needs this large pin/pout difference

	PLP lp;
	Partition zeta = lp.run(G);

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	DEBUG("modularity produced by LabelPropagation: " , mod);
	DEBUG("number of clusters produced by LabelPropagation: k=" , zeta.numberOfSubsets());

	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta)) << "the resulting partition should be a proper clustering";
	EXPECT_TRUE(GraphClusteringTools::equalClusterings(zeta, reference, G)) << "LP should detect exactly the reference clustering";

}



TEST_F(ClusteringAlgoGTest, testLabelPropagationOnDisconnectedGraph) {
	GraphGenerator graphGenerator;
	int n = 100;
	int k = 2; // number of clusters
	Graph G = graphGenerator.makeClusteredRandomGraph(n, k, 1.0, 0.0);

	PLP lp;
	Partition zeta = lp.run(G);

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	DEBUG("modularity produced by LabelPropagation: " , mod);

	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta)) << "the resulting partition should be a proper clustering";
	EXPECT_EQ(k, zeta.numberOfSubsets()) << " " << k << " clusters are easy to detect"; //FIXME

}


TEST_F(ClusteringAlgoGTest, testLabelPropagationOnSingleNodeWithSelfLoop) {
	Graph G(1);
	node v = 0;
	G.setWeight(v, v, 42.0);

	PLP lp;
	Partition zeta = lp.run(G);

	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta));
	EXPECT_TRUE(GraphClusteringTools::isSingletonClustering(G, zeta));
	EXPECT_TRUE(GraphClusteringTools::isOneClustering(G, zeta)); //FIXME does this make sense? singleton and one partition at the same time.

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	DEBUG("modularity produced by LabelPropagation: " , mod);

}



TEST_F(ClusteringAlgoGTest, testLabelPropagationOnManySmallClusters) {
	int64_t n = 1000;
	int k = 100; // number of clusters
	double pin = 1.0;
	double pout = 0.0;

	GraphGenerator graphGen;
	std::pair<Graph, Partition> G_ref = graphGen.makeClusteredRandomGraphWithReferenceClustering(n, k, pin, pout);


	PLP lp;
	Partition zeta = lp.run(G_ref.first);

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G_ref.first);
	DEBUG("modularity produced by LabelPropagation: " , mod);
	DEBUG("number of clusters produced by LabelPropagation: k=" , zeta.numberOfSubsets());

	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G_ref.first, zeta)) << "the resulting partition should be a proper clustering";
	EXPECT_TRUE(GraphClusteringTools::equalClusterings(zeta, G_ref.second, G_ref.first)) << "Can LabelPropagation detect the reference clustering?";

}

TEST_F(ClusteringAlgoGTest, testLouvain) {
	count n = 500;
	count k = 25;
	double pin = 0.9;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	PLM louvain;
	Partition zeta = louvain.run(G);

	INFO("number of clusters: " , zeta.numberOfSubsets());

	Modularity modularity;
	INFO("modularity: " , modularity.getQuality(zeta, G));

}


TEST_F(ClusteringAlgoGTest, testLouvainParallelSimple) {
	count n = 500;
	count k = 25;
	double pin = 0.9;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	PLM louvain("simple");
	Partition zeta = louvain.run(G);

	INFO("number of clusters: " , zeta.numberOfSubsets());

	Modularity modularity;
	INFO("modularity: " , modularity.getQuality(zeta, G));

}

/*
TEST_F(ClusteringAlgoGTest, testLouvainParallel2Naive) {
	count n = 1000;
	count k = 100;
	double pin = 1.0;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	LouvainParallel louvain;
	Partition zeta = louvain.run(G);

	INFO("number of clusters: " , zeta.numberOfSubsets());

	Modularity modularity;
	INFO("modularity: " , modularity.getQuality(zeta, G));

}
*/


TEST_F(ClusteringAlgoGTest, testLouvainParallelBalanced) {
	count n = 500;
	count k = 25;
	double pin = 0.9;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	PLM louvain("balanced");
	Partition zeta = louvain.run(G);

	INFO("number of clusters: " , zeta.numberOfSubsets());

	Modularity modularity;
	INFO("modularity: " , modularity.getQuality(zeta, G));

}


TEST_F(ClusteringAlgoGTest, testCNM) {
	count n = 200;
	count k = 25;
	double pin = 0.9;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	Modularity modularity;

	// CNM with PQ
	CNM cnm;
	Partition clustering = cnm.run(G);
	INFO("CNM number of clusters: " , clustering.numberOfSubsets());
	INFO("modularity clustered random graph: " , modularity.getQuality(clustering, G));
	// EXPECT_GE(modularity.getQuality(clustering, G), 0.5);
	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, clustering));

}


TEST_F(ClusteringAlgoGTest, testCNMandLouvain) {
	Modularity modularity;
	CNM cnm;
	PLM louvain;
	METISGraphReader reader;
	Graph jazz = reader.read("input/jazz.graph");
	// this takes much longer than a unit test should
	// Graph blog = reader.read("input/polblogs.graph");



	// *** jazz graph
	// Louvain
	Partition clustering = louvain.run(jazz);
	INFO("Louvain number of jazz clusters: " , clustering.numberOfSubsets());
	INFO("Louvain modularity jazz graph: " , modularity.getQuality(clustering, jazz));

	// CNM
	clustering = cnm.run(jazz);
	INFO("CNM number of jazz clusters: " , clustering.numberOfSubsets());
	INFO("CNM modularity jazz graph: " , modularity.getQuality(clustering, jazz));

//	// *** blog graph
//	// CNM
//	clustering = cnm.run(blog);
//	INFO("CNM number of blog clusters: " , clustering.numberOfSubsets());
//	INFO("CNM modularity blog graph: " , modularity.getQuality(clustering, jazz));
//
//	// Louvain
//	clustering = louvain.run(blog);
//	INFO("Louvain number of blog clusters: " , clustering.numberOfSubsets());
//	INFO("Louvain modularity blog graph: " , modularity.getQuality(clustering, jazz));
}


TEST_F(ClusteringAlgoGTest, testParallelAgglomerativeAndPLM2) {
	Modularity modularity;
	ParallelAgglomerativeClusterer aggl;
	PLM2 louvain;
	METISGraphReader reader;
	Graph jazz = reader.read("input/jazz.graph");
	Graph blog = reader.read("input/polblogs.graph");


	// *** jazz graph
	// parallel agglomerative
	Partition clustering = aggl.run(jazz);
	INFO("Match-AGGL number of jazz clusters: " , clustering.numberOfSubsets());
	INFO("Match-AGGL modularity jazz graph:   " , modularity.getQuality(clustering, jazz));

	// Louvain
	clustering = louvain.run(jazz);
	INFO("Louvain number of jazz clusters: " , clustering.numberOfSubsets());
	INFO("Louvain modularity jazz graph:   " , modularity.getQuality(clustering, jazz));


	// *** blog graph
	// parallel agglomerative
	clustering = aggl.run(blog);
	INFO("Match-AGGL number of blog clusters: " , clustering.numberOfSubsets());
	INFO("Match-AGGL modularity blog graph:   " , modularity.getQuality(clustering, blog));

	// Louvain
	clustering = louvain.run(blog);
	INFO("Louvain number of blog clusters: " , clustering.numberOfSubsets());
	INFO("Louvain modularity blog graph:   " , modularity.getQuality(clustering, blog));
}




TEST_F(ClusteringAlgoGTest, testEPPFactory) {

	EPPFactory factory;
	EPP epp = factory.make(4, "PLP", "PLM");

	METISGraphReader reader;
	Graph jazz = reader.read("input/jazz.graph");
	Partition zeta = epp.run(jazz);

	INFO("number of clusters: " , zeta.numberOfSubsets());

	EXPECT_TRUE(GraphClusteringTools::isProperClustering(jazz, zeta));
}



TEST_F(ClusteringAlgoGTest, testPLM2) {
	METISGraphReader reader;
	Modularity modularity;
	Graph G = reader.read("input/PGPgiantcompo.graph");

	PLM2 plm(false, 1.0);
	Partition zeta = plm.run(G);

	INFO("number of clusters: " , zeta.numberOfSubsets());
	INFO("modularity: " , modularity.getQuality(zeta, G));
	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta));

	PLM2 plmr(true, 1.0);
	Partition zeta2 = plmr.run(G);

	INFO("number of clusters: " , zeta2.numberOfSubsets());
	INFO("modularity: " , modularity.getQuality(zeta2, G));
	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta2));

}

TEST_F(ClusteringAlgoGTest, testCommunityGraph) {
	CommunityGraph com;
	METISGraphReader reader;
	Graph G = reader.read("input/jazz.graph");
	ClusteringGenerator clusteringGen;

	Partition one = clusteringGen.makeOneClustering(G);
	com.run(G, one);
	EXPECT_EQ(1, com.getGraph().numberOfNodes());

	Partition singleton = clusteringGen.makeSingletonClustering(G);
	com.run(G, singleton);
	EXPECT_EQ(G.numberOfNodes(), com.getGraph().numberOfNodes());

	Partition zeta = (new PLP())->run(G);
	com.run(G, zeta);
	EXPECT_EQ(zeta.numberOfSubsets(), com.getGraph().numberOfNodes());
}


} /* namespace NetworKit */

#endif /*NOGTEST */
