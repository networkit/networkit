/*
 * ClusteringTest.cpp
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "ClusteringGTest.h"

#include "../../auxiliary/Log.h"
#include "../Clustering.h"
#include "../Modularity.h"
#include "../ModularitySequential.h"
#include "../Coverage.h"
#include "../ClusteringGenerator.h"
#include "../JaccardMeasure.h"
#include "../NodeStructuralRandMeasure.h"
#include "../GraphStructuralRandMeasure.h"
#include "../../graph/GraphGenerator.h"
//#include "../../io/METISGraphReader.h"
//#include "../../io/ClusteringReader.h"
#include "../NMIDistance.h"
#include "../DynamicNMIDistance.h"
#include "../../auxiliary/NumericTools.h"
#include "../../dcd/DynCDSetup.h"
#include "../../generators/DynamicBarabasiAlbertGenerator.h"
//#include "../../dcd/DynamicPLP.h"
#include "../../dcd/TDynamicLabelPropagation.h"

namespace NetworKit {



TEST_F(ClusteringGTest, testModularity) {
	GraphGenerator graphGenerator;

	count n = 100;

	DEBUG("testing modularity on clustering of complete graph with " , n , " nodes");


	Graph G = graphGenerator.makeCompleteGraph(n);
	DEBUG("total edge weight: " , G.totalEdgeWeight());

	ClusteringGenerator clusteringGenerator;

	Clustering singleton = clusteringGenerator.makeSingletonClustering(G);
	Clustering one = clusteringGenerator.makeOneClustering(G);

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

TEST_F(ClusteringGTest, testCoverage) {
	GraphGenerator graphGenerator;

	count n = 100;

	DEBUG("testing coverage on clustering of complete graph with " , n , " nodes");


	Graph G = graphGenerator.makeCompleteGraph(n);

	ClusteringGenerator clusteringGenerator;

	Clustering singleton = clusteringGenerator.makeSingletonClustering(G);
	Clustering one = clusteringGenerator.makeOneClustering(G);

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


TEST_F(ClusteringGTest, testClusteringEquality) {
	count n = 100;
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
	count n = 100;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Clustering singleton = clusteringGen.makeSingletonClustering(G);
	Clustering random = clusteringGen.makeRandomClustering(G, 10);

	JaccardMeasure jaccard;
	double j = jaccard.getDissimilarity(G, singleton, random);

	EXPECT_EQ(1.0, j) << "The singleton clustering and any other clustering compare with a dissimilarity of 1.0";

}


TEST_F(ClusteringGTest, testNodeStructuralRandMeasure) {
	count n = 100;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Clustering one1 = clusteringGen.makeOneClustering(G);
	Clustering one2 = clusteringGen.makeOneClustering(G);

	NodeStructuralRandMeasure rand;
	double r = rand.getDissimilarity(G, one1, one2);

	EXPECT_EQ(0.0, r) << "Identical clusterings should compare with a dissimilarity of 0.0";

}

TEST_F(ClusteringGTest, testGraphStructuralRandMeasure) {
	count n = 100;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Clustering one1 = clusteringGen.makeOneClustering(G);
	Clustering one2 = clusteringGen.makeOneClustering(G);

	GraphStructuralRandMeasure rand;
	double r = rand.getDissimilarity(G, one1, one2);

	EXPECT_EQ(0.0, r) << "Identical clusterings should compare with a dissimilarity of 0.0";

}


TEST_F(ClusteringGTest, testCompact) {
	count n = 50;
	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	ClusteringGenerator clusteringGen;
	Clustering clustering = clusteringGen.makeRandomClustering(G, 2*n);
	clustering.compact();

	EXPECT_LE(clustering.upperBound(), n);
}

TEST_F(ClusteringGTest, testModularityParallelVsSequential) {

	Modularity modularityPar;
	ModularitySequential modularitySeq;

	count n = 500;
	GraphGenerator graphGen;
	ClusteringGenerator clusteringGen;

	INFO("making random graph");
	Graph G = graphGen.makeRandomGraph(n, 0.2);
	INFO("making random clustering");
	Clustering zeta = clusteringGen.makeRandomClustering(G, 42);

	INFO("calculating modularity in parallel");
	double modPar = modularityPar.getQuality(zeta, G);
	INFO("calculating modularity sequentially");
	double modSeq = modularitySeq.getQuality(zeta, G);

	// EXPECT_EQ(modPar, modSeq) << "Modularity values should be equal no matter if calculated in parallel or sequentially";
	EXPECT_TRUE(Aux::NumericTools::equal(modPar, modSeq, 1e-12)) << "Modularity values should be equal within a small error no matter if calculated in parallel or sequentially";

}



//TEST_F(ClusteringGTest, testModularityParallelVsSequentialOnLargeGraph) {
//
//	Modularity modularityPar;
//	ModularitySequential modularitySeq;
//
//	ClusteringGenerator clusteringGen;
//
//	METISGraphReader reader;
//	Graph G = reader.read("graphs/Benchmark/uk-2002.graph"); // FIXME: hardcoded file name
//	Clustering zeta = clusteringGen.makeRandomClustering(G, 42);
//
//	double modPar = modularityPar.getQuality(zeta, G);
//	double modSeq = modularitySeq.getQuality(zeta, G);
//
//	EXPECT_EQ(modPar, modSeq) << "Modularity values should be equal no matter if calculated in parallel or sequentially";
//
//}


//TEST_F(ClusteringGTest, testModularityWithStoredClustering) {
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
//	Clustering zeta = clusteringReader.read(clusteringPath);
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



TEST_F(ClusteringGTest, testNMIDistance) {
	// two 1-clusterings should have NMID = 0 because H is 0
	GraphGenerator gen;
	Graph G = gen.makeErdosRenyiGraph(10, 1.0);

	ClusteringGenerator clustGen;
	Clustering one1 = clustGen.makeOneClustering(G);
	Clustering one2 = clustGen.makeOneClustering(G);

	NMIDistance NMID;
	double distOne = NMID.getDissimilarity(G, one1, one2);

	INFO("NMID for two 1-clusterings: " , distOne);
	EXPECT_TRUE(Aux::NumericTools::equal(0.0, distOne)) << "NMID of two 1-clusterings should be 0.0";


	Clustering singleton1 = clustGen.makeSingletonClustering(G);
	Clustering singleton2 = clustGen.makeSingletonClustering(G);

	double distSingleton = NMID.getDissimilarity(G, singleton1, singleton2);
	INFO("NMID for two singleton clusterings: " , distSingleton);


	EXPECT_TRUE(Aux::NumericTools::equal(0.0, distSingleton)) << "NMID of two identical singleton clusterings should be 0.0";

	Clustering random1 = clustGen.makeRandomClustering(G, 2);
	Clustering random2 = clustGen.makeRandomClustering(G, 2);

	double distRandom = NMID.getDissimilarity(G, random1, random2);
	INFO("NMID for two random clusterings: " , distRandom);

}

TEST_F(ClusteringGTest, tryDynamicNMIDistance) {
	// FIXME: sometimes throws an assertion error
	// two 1-clusterings should have NMID = 0 because H is 0
	GraphGenerator gen;
	Graph G = gen.makeErdosRenyiGraph(10, 1.0);

	ClusteringGenerator clustGen;
	Clustering one1 = clustGen.makeOneClustering(G);
	Clustering one2 = clustGen.makeOneClustering(G);

	DynamicNMIDistance dynNMID;
	double distOne = dynNMID.getDissimilarity(G, one1, one2);

	INFO("Dyn NMID for two 1-clusterings: " , distOne);
	EXPECT_TRUE(Aux::NumericTools::equal(distOne, 0.0)) << "Dyn NMID of two 1-clusterings should be 0.0";


	Clustering singleton1 = clustGen.makeSingletonClustering(G);
	Clustering singleton2 = clustGen.makeSingletonClustering(G);

	double distSingleton = dynNMID.getDissimilarity(G, singleton1, singleton2);
	INFO("Dyn NMID for two singleton clusterings: " , distSingleton);

	EXPECT_TRUE(Aux::NumericTools::equal(distSingleton, 0.0)) << "Dyn NMID of two identical singleton clusterings should be 0.0";

	Clustering random1 = clustGen.makeRandomClustering(G, 2);
	Clustering random2 = clustGen.makeRandomClustering(G, 2);

	double distRandom = dynNMID.getDissimilarity(G, random1, random2);
	INFO("Dyn NMID for two random clusterings: " , distRandom);


	// now dynamic graph(s)
	DynamicGraphSource* dynGen = new DynamicBarabasiAlbertGenerator(1);
	DynamicCommunityDetector* dynLP1 = new TDynamicLabelPropagation<Isolate>();

	std::vector<DynamicCommunityDetector*> detectors = { dynLP1 };
	count deltaT = 6;
	count tMax = 10 * deltaT;
	DynCDSetup setup(*dynGen, detectors, tMax, deltaT);

	setup.run();

	G = setup.getGraphCopy();
	std::vector<Clustering>& myresults = setup.dynamicClusteringTimelines[0];
	Clustering& currentClustering = myresults.back();
	Clustering& oldClustering = myresults[0];
	EXPECT_TRUE(currentClustering.isProper(G)) << "clustering in the sequence should be a proper clustering of G";

	double distSeq = dynNMID.getDissimilarity(G, oldClustering, currentClustering);
	INFO("Dyn NMID for last and first clustering: " , distSeq);

}


TEST_F(ClusteringGTest, testNumberOfClusters) {
	GraphGenerator graphGenerator;
	count n = 42;
	count k = 4; // number of clusters
	Graph G = graphGenerator.makeCompleteGraph(n);

	ClusteringGenerator clusteringGenerator;
	Clustering singleton = clusteringGenerator.makeSingletonClustering(G);
	Clustering one = clusteringGenerator.makeOneClustering(G);

	EXPECT_EQ(n, singleton.numberOfClusters());
	EXPECT_EQ(1, one.numberOfClusters());
}


TEST_F(ClusteringGTest, testGetMembers) {
	GraphGenerator graphGenerator;
	count n = 42;
	Graph G = graphGenerator.makeCompleteGraph(n);
	ClusteringGenerator clusteringGenerator;
	Clustering one = clusteringGenerator.makeOneClustering(G);
	node v = 0; // arbitrary node
	cluster C = one[v]; // get cluster id
	std::vector<node> members = one.getMembers(C);

	EXPECT_EQ(n, members.size()) << "all nodes must be in the 1-cluster";

}


TEST_F(ClusteringGTest, testClusterSizes) {
	GraphGenerator graphGenerator;
	count n = 42;
	Graph G = graphGenerator.makeCompleteGraph(n);
	ClusteringGenerator clusteringGenerator;
	Clustering one = clusteringGenerator.makeOneClustering(G);

	std::vector<count> sizes = one.clusterSizes();

	EXPECT_EQ(1, sizes.size()) << "only one entry for 1-clustering";
	EXPECT_EQ(n, sizes.at(0)) << "size of only cluster should be n";
}




} /* namespace NetworKit */

#endif /*NOGTEST */
