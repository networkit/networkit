/*
 * CommunityGTest.cpp
 *
 *  Created on: 27.02.2014
 *      Author: cls
 */

#include "CommunityGTest.h"

#include "../PLP.h"
#include "../PLM.h"
#include "../CNM.h"
#include "../ParallelAgglomerativeClusterer.h"
#include "../../community/Modularity.h"
#include "../../graph/GraphGenerator.h"
#include "../../community/ClusteringGenerator.h"
#include "../../io/METISGraphReader.h"
#include "../EPP.h"
#include "../../overlap/HashingOverlapper.h"
#include "../EPPFactory.h"
#include "../CommunityGraph.h"
#include "../PLM.h"
#include "../PLMOld.h"
#include "../../community/GraphClusteringTools.h"
#include "../../auxiliary/Log.h"
#include "../../structures/Partition.h"
#include "../Modularity.h"
#include "../ModularitySequential.h"
#include "../Coverage.h"
#include "../ClusteringGenerator.h"
#include "../JaccardMeasure.h"
#include "../NodeStructuralRandMeasure.h"
#include "../GraphStructuralRandMeasure.h"
#include "../../graph/GraphGenerator.h"
#include "../NMIDistance.h"
#include "../DynamicNMIDistance.h"
#include "../../auxiliary/NumericTools.h"
#include "../../generators/DynamicBarabasiAlbertGenerator.h"
#include "../SampledGraphStructuralRandMeasure.h"
#include "../SampledNodeStructuralRandMeasure.h"
#include "../../community/GraphClusteringTools.h"

namespace NetworKit {

TEST_F(CommunityGTest, testEnsemblePreprocessing) {
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



TEST_F(CommunityGTest, testLabelPropagationOnUniformGraph) {
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


TEST_F(CommunityGTest, testLabelPropagationOnClusteredGraph_ForNumberOfClusters) {
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


TEST_F(CommunityGTest, testLabelPropagationOnClusteredGraph_ForEquality) {
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



TEST_F(CommunityGTest, testLabelPropagationOnDisconnectedGraph) {
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


TEST_F(CommunityGTest, testLabelPropagationOnSingleNodeWithSelfLoop) {
	Graph G(1, true);
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



TEST_F(CommunityGTest, testLabelPropagationOnManySmallClusters) {
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

TEST_F(CommunityGTest, testLouvain) {
	count n = 500;
	count k = 25;
	double pin = 0.9;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	PLMOld louvain;
	Partition zeta = louvain.run(G);

	INFO("number of clusters: " , zeta.numberOfSubsets());

	Modularity modularity;
	INFO("modularity: " , modularity.getQuality(zeta, G));

}


TEST_F(CommunityGTest, testLouvainParallelSimple) {
	count n = 500;
	count k = 25;
	double pin = 0.9;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	PLMOld louvain("simple");
	Partition zeta = louvain.run(G);

	INFO("number of clusters: " , zeta.numberOfSubsets());

	Modularity modularity;
	INFO("modularity: " , modularity.getQuality(zeta, G));

}

/*
TEST_F(CommunityGTest, testLouvainParallel2Naive) {
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


TEST_F(CommunityGTest, testLouvainParallelBalanced) {
	count n = 500;
	count k = 25;
	double pin = 0.9;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	PLMOld louvain("balanced");
	Partition zeta = louvain.run(G);

	INFO("number of clusters: " , zeta.numberOfSubsets());

	Modularity modularity;
	INFO("modularity: " , modularity.getQuality(zeta, G));

}


TEST_F(CommunityGTest, testCNMandLouvainRandom) {
	count n = 400;
	count k = 20;
	double pin = 0.9;
	double pout = 0.0005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	Modularity modularity;

	// CNM with PQ
	CNM cnm;
	Partition clustering = cnm.run(G);
	INFO("CNM number of clusters: " , clustering.numberOfSubsets());
	INFO("modularity clustered random graph: " , modularity.getQuality(clustering, G));
	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, clustering));

	// Louvain
	PLM louvain;
	clustering = louvain.run(G);
	INFO("Louvain number of clusters: " , clustering.numberOfSubsets());
	INFO("modularity clustered random graph: " , modularity.getQuality(clustering, G));
	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, clustering));
}


TEST_F(CommunityGTest, testCNMandLouvainReal) {
	Modularity modularity;
	CNM cnm;
	PLM louvain;
	METISGraphReader reader;
	Graph jazz = reader.read("input/jazz.graph");

	// CNM
	Partition clustering = cnm.run(jazz);
	INFO("CNM number of jazz clusters: " , clustering.numberOfSubsets());
	INFO("CNM modularity jazz graph: " , modularity.getQuality(clustering, jazz));

	// Louvain
	clustering = louvain.run(jazz);
	INFO("Louvain number of jazz clusters: " , clustering.numberOfSubsets());
	INFO("Louvain modularity jazz graph: " , modularity.getQuality(clustering, jazz));
}


TEST_F(CommunityGTest, testParallelAgglomerativeAndLouvain) {
	Modularity modularity;
	ParallelAgglomerativeClusterer aggl;
	PLMOld louvain;
	METISGraphReader reader;
	Graph jazz = reader.read("input/jazz.graph");
	Graph blog = reader.read("input/polblogs.graph");


	// *** jazz graph
	// aggl
	Partition clustering = aggl.run(jazz);
	INFO("Match-AGGL number of jazz clusters: " , clustering.numberOfSubsets());
	INFO("Match-AGGL modularity jazz graph:   " , modularity.getQuality(clustering, jazz));

	// Louvain
	clustering = louvain.run(jazz);
	INFO("Louvain number of jazz clusters: " , clustering.numberOfSubsets());
	INFO("Louvain modularity jazz graph:   " , modularity.getQuality(clustering, jazz));


	// *** blog graph
	// CNM
	clustering = aggl.run(blog);
	INFO("Match-AGGL number of blog clusters: " , clustering.numberOfSubsets());
	INFO("Match-AGGL modularity blog graph:   " , modularity.getQuality(clustering, blog));

	// Louvain
	clustering = louvain.run(blog);
	INFO("Louvain number of blog clusters: " , clustering.numberOfSubsets());
	INFO("Louvain modularity blog graph:   " , modularity.getQuality(clustering, blog));
}




TEST_F(CommunityGTest, testEPPFactory) {

	EPPFactory factory;
	EPP epp = factory.make(4, "PLP", "PLM");

	METISGraphReader reader;
	Graph jazz = reader.read("input/jazz.graph");
	Partition zeta = epp.run(jazz);

	INFO("number of clusters: " , zeta.numberOfSubsets());

	EXPECT_TRUE(GraphClusteringTools::isProperClustering(jazz, zeta));
}



TEST_F(CommunityGTest, testPLM) {
	METISGraphReader reader;
	Modularity modularity;
	Graph G = reader.read("input/PGPgiantcompo.graph");

	PLM plm(false, 1.0);
	Partition zeta = plm.run(G);

	INFO("number of clusters: " , zeta.numberOfSubsets());
	INFO("modularity: " , modularity.getQuality(zeta, G));
	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta));

	PLM plmr(true, 1.0);
	Partition zeta2 = plmr.run(G);

	INFO("number of clusters: " , zeta2.numberOfSubsets());
	INFO("modularity: " , modularity.getQuality(zeta2, G));
	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, zeta2));

}

TEST_F(CommunityGTest, testCommunityGraph) {
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


TEST_F(CommunityGTest, testModularity) {
	GraphGenerator graphGenerator;

	count n = 100;

	DEBUG("testing modularity on clustering of complete graph with " , n , " nodes");


	Graph G = graphGenerator.makeCompleteGraph(n);
	DEBUG("total edge weight: " , G.totalEdgeWeight());

	ClusteringGenerator clusteringGenerator;

	Partition singleton = clusteringGenerator.makeSingletonClustering(G);
	Partition one = clusteringGenerator.makeOneClustering(G);

	Modularity modularity;

	DEBUG("calculating modularity for singleton clustering");
	double modSingleton = modularity.getQuality(singleton, G);

	DEBUG("calculating modularity for 1-clustering");
	double modOne = modularity.getQuality(one, G);

	DEBUG("mod(singleton-clustering) = " , modSingleton);
	DEBUG("mod(1-clustering) = " , modOne);


	EXPECT_EQ(0.0, modOne) << "1-clustering should have modularity of 0.0";
	EXPECT_GE(0.0, modSingleton) << "singleton clustering should have modularity less than 0.0";

}

TEST_F(CommunityGTest, testCoverage) {
	GraphGenerator graphGenerator;

	count n = 100;

	DEBUG("testing coverage on clustering of complete graph with " , n , " nodes");


	Graph G = graphGenerator.makeCompleteGraph(n);

	ClusteringGenerator clusteringGenerator;

	Partition singleton = clusteringGenerator.makeSingletonClustering(G);
	Partition one = clusteringGenerator.makeOneClustering(G);

	Coverage coverage;

	DEBUG("calculating coverage for singleton clustering");
	double covSingleton = coverage.getQuality(singleton, G);

	DEBUG("calculating coverage for 1-clustering");
	double covOne = coverage.getQuality(one, G);

	DEBUG("mod(singleton-clustering) = " , covSingleton);
	DEBUG("mod(1-clustering) = " , covOne);


	EXPECT_EQ(1.0, covOne) << "1-clustering should have coverage of 1.0";
	EXPECT_GE(0.0, covSingleton) << "singleton clustering should have coverage of 0.0";

}

// TODO necessary testcase? move equals to some class ?
TEST_F(CommunityGTest, testClusteringEquality) {
	count n = 100;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Partition one1 = clusteringGen.makeOneClustering(G);
	Partition one2 = clusteringGen.makeOneClustering(G);

	EXPECT_TRUE(GraphClusteringTools::equalClusterings(one1, one2, G)) << "two 1-clusterings of G should be equal";

	Partition singleton1 = clusteringGen.makeSingletonClustering(G);
	Partition singleton2 = clusteringGen.makeSingletonClustering(G);

	EXPECT_TRUE(GraphClusteringTools::equalClusterings(singleton1, singleton2, G)) << "two singleton clusterings of G should be equal";

}



TEST_F(CommunityGTest, testJaccardMeasure) {
	count n = 100;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Partition singleton = clusteringGen.makeSingletonClustering(G);
	Partition random = clusteringGen.makeRandomClustering(G, 10);

	JaccardMeasure jaccard;
	double j = jaccard.getDissimilarity(G, singleton, random);

	EXPECT_EQ(1.0, j) << "The singleton clustering and any other clustering compare with a dissimilarity of 1.0";

}


TEST_F(CommunityGTest, testNodeStructuralRandMeasure) {
	count n = 100;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Partition one1 = clusteringGen.makeOneClustering(G);
	Partition one2 = clusteringGen.makeOneClustering(G);

	NodeStructuralRandMeasure rand;
	double r = rand.getDissimilarity(G, one1, one2);

	EXPECT_EQ(0.0, r) << "Identical clusterings should compare with a dissimilarity of 0.0";

}

TEST_F(CommunityGTest, testGraphStructuralRandMeasure) {
	count n = 100;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Partition one1 = clusteringGen.makeOneClustering(G);
	Partition one2 = clusteringGen.makeOneClustering(G);

	GraphStructuralRandMeasure rand;
	double r = rand.getDissimilarity(G, one1, one2);

	EXPECT_EQ(0.0, r) << "Identical clusterings should compare with a dissimilarity of 0.0";

}


TEST_F(CommunityGTest, testModularityParallelVsSequential) {

	Modularity modularityPar;
	ModularitySequential modularitySeq;

	count n = 500;
	GraphGenerator graphGen;
	ClusteringGenerator clusteringGen;

	INFO("making random graph");
	Graph G = graphGen.makeRandomGraph(n, 0.2);
	INFO("making random clustering");
	Partition zeta = clusteringGen.makeRandomClustering(G, 42);

	INFO("calculating modularity in parallel");
	double modPar = modularityPar.getQuality(zeta, G);
	INFO("calculating modularity sequentially");
	double modSeq = modularitySeq.getQuality(zeta, G);

	// EXPECT_EQ(modPar, modSeq) << "Modularity values should be equal no matter if calculated in parallel or sequentially";
	EXPECT_TRUE(Aux::NumericTools::equal(modPar, modSeq, 1e-12)) << "Modularity values should be equal within a small error no matter if calculated in parallel or sequentially";

}



//TEST_F(CommunityGTest, testModularityParallelVsSequentialOnLargeGraph) {
//
//	Modularity modularityPar;
//	ModularitySequential modularitySeq;
//
//	ClusteringGenerator clusteringGen;
//
//	METISGraphReader reader;
//	Graph G = reader.read("graphs/Benchmark/uk-2002.graph"); // FIXME: hardcoded file name
//	Partition zeta = clusteringGen.makeRandomClustering(G, 42);
//
//	double modPar = modularityPar.getQuality(zeta, G);
//	double modSeq = modularitySeq.getQuality(zeta, G);
//
//	EXPECT_EQ(modPar, modSeq) << "Modularity values should be equal no matter if calculated in parallel or sequentially";
//
//}


//TEST_F(CommunityGTest, testModularityWithStoredClustering) {
//
//	std::string graphPath;
//	std::cout << "[INPUT] .graph file path >" << std::endl;
//	std::getline(std::cin, graphPath);
//
//	std::string clusteringPath;
//	std::cout << "[INPUT] .clust/.ptn file path >" << std::endl;
//	std::getline(std::cin, clusteringPath);
//
//	std::string evalPath;
//	std::cout << "[INPUT] .eval file path >" << std::endl;
//	std::getline(std::cin, evalPath);
//
//	INFO("reading graph from: " , graphPath);
//	METISGraphReader graphReader;
//	Graph G = graphReader.read(graphPath);
//
//	ClusteringReader clusteringReader;
//	INFO("reading clustering from: " , clusteringPath);
//	Partition zeta = clusteringReader.read(clusteringPath);
//
//	INFO("reading modularity value from .eval file: " , evalPath);
//	std::ifstream evalFile(evalPath);
//	std::string evalLine;
//	std::getline(evalFile, evalLine);
//	double evalMod = std::atof(evalLine.c_str());
//	INFO("modularity from .eval file: " , evalMod);
//
//	Modularity modularity;
//	INFO("calculating modularity in parallel");
//	double modPar = modularity.getQuality(zeta, G);
//	INFO("modPar: " , modPar);
//
//	ModularitySequential modularitySeq;
//	INFO("calculating modularity sequentially");
//	double modSeq = modularitySeq.getQuality(zeta, G);
//	INFO("modSeq: " , modSeq);
//
//	EXPECT_EQ(modSeq, modPar) << "Modularity values should be equal no matter if calculated in parallel or sequentially";
//	EXPECT_EQ(modSeq, evalMod) << "modSeq should be agree with DIMACS challenge evaluation";
//	EXPECT_EQ(modPar, evalMod) << "modPar should be agree with DIMACS challenge evaluation";
//
//}



TEST_F(CommunityGTest, testNMIDistance) {
	// two 1-clusterings should have NMID = 0 because H is 0
	GraphGenerator gen;
	Graph G = gen.makeErdosRenyiGraph(10, 1.0);

	ClusteringGenerator clustGen;
	Partition one1 = clustGen.makeOneClustering(G);
	Partition one2 = clustGen.makeOneClustering(G);

	NMIDistance NMID;
	double distOne = NMID.getDissimilarity(G, one1, one2);

	INFO("NMID for two 1-clusterings: " , distOne);
	EXPECT_TRUE(Aux::NumericTools::equal(0.0, distOne)) << "NMID of two 1-clusterings should be 0.0";


	Partition singleton1 = clustGen.makeSingletonClustering(G);
	Partition singleton2 = clustGen.makeSingletonClustering(G);

	double distSingleton = NMID.getDissimilarity(G, singleton1, singleton2);
	INFO("NMID for two singleton clusterings: " , distSingleton);


	EXPECT_TRUE(Aux::NumericTools::equal(0.0, distSingleton)) << "NMID of two identical singleton clusterings should be 0.0";

	Partition random1 = clustGen.makeRandomClustering(G, 2);
	Partition random2 = clustGen.makeRandomClustering(G, 2);

	double distRandom = NMID.getDissimilarity(G, random1, random2);
	INFO("NMID for two random clusterings: " , distRandom);

}


TEST_F(CommunityGTest, testSampledRandMeasures) {
	GraphGenerator graphGenerator;
	count n = 42;
	Graph G = graphGenerator.makeCompleteGraph(n);
	ClusteringGenerator clusteringGenerator;
	Partition one = clusteringGenerator.makeOneClustering(G);
	Partition singleton = clusteringGenerator.makeSingletonClustering(G);

	SampledNodeStructuralRandMeasure nRand(20);
	SampledGraphStructuralRandMeasure gRand(20);

	double nDis = nRand.getDissimilarity(G, one, singleton);
	double gDis = gRand.getDissimilarity(G, one, singleton);

	DEBUG("node structural dissimilarity: ", nDis);
	DEBUG("graph structural dissimilarity: ", gDis);
}


TEST_F(CommunityGTest, testParallelAgglomerativeAndPLM) {
	Modularity modularity;
	ParallelAgglomerativeClusterer aggl;
	PLM louvain;
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



} /* namespace NetworKit */
