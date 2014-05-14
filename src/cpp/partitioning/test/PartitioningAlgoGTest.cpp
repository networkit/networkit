/*
 * PartitioningAlgoGTest.cpp
 *
 *  Created on: 19.06.2013
 *      Author: Henning Meyerhenke (meyerhenke@kit.edu)
 */

#include "PartitioningAlgoGTest.h"
#include "../../community/GraphClusteringTools.h"

#ifndef NOGTEST

namespace NetworKit {

TEST_F(PartitioningAlgoGTest, tryBalancedLabelPropagationOnClusteredGraph) {
	GraphGenerator graphGenerator;
	int64_t n = 300;
	count k = 3; // number of clusters
	Graph G = graphGenerator.makeClusteredRandomGraph(n, k, 1.0, 0.001);

	double exponent = 1.75;
	BalancedLabelPropagation lp(exponent);
	ClusteringGenerator gen;
	Partition zeta = gen.makeRandomClustering(G, k);
	zeta = lp.rerun(G, k, zeta);

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	EdgeCut edgeCut;
	double cut = edgeCut.getQuality(zeta, G);

	DEBUG("modularity produced by BalancedLabelPropagation: " , mod , ", cut: " , cut);

	EXPECT_TRUE(GraphClusteringTools::isProperClustering(G,zeta)) << "the resulting partition should be a proper clustering";
	EXPECT_EQ(k, zeta.numberOfSubsets()) << " " << k << " clusters are easy to detect";

	std::vector<count> clusterSizes = zeta.subsetSizes();
	for (index p = 0; p < k; ++p) {
		DEBUG("size of cluster " , p , ": " , clusterSizes[p]);
	}
}


TEST_F(PartitioningAlgoGTest, tryBalancedLabelPropagationOnRealGraph) {
	Modularity modularity;
	EdgeCut edgeCut;
	METISGraphReader reader;
	Graph airfoil1 = reader.read("input/airfoil1.graph");

	double exponent = 1.75;
	BalancedLabelPropagation blp(exponent);
	count k = 4;

	// *** airfoil1 graph
	// BLP
	ClusteringGenerator gen;
	Partition partition = gen.makeRandomClustering(airfoil1, k);
	partition = blp.rerun(airfoil1, k, partition);
	double cut = edgeCut.getQuality(partition, airfoil1);
	INFO("BLP number of airfoil1 clusters: " , partition.numberOfSubsets());
	INFO("BLP modularity airfoil1 graph:   " , modularity.getQuality(partition, airfoil1) , ", cut: " , cut);
	std::vector<count> clusterSizes = partition.subsetSizes();
	for (index p = 0; p < k; ++p) {
		DEBUG("size of cluster " , p , ": " , clusterSizes[p]);
	}
}

TEST_F(PartitioningAlgoGTest, tryMultilevelBalancedLabelPropagationOnRealGraph) {
	Modularity modularity;
	EdgeCut edgeCut;
#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
	DibapGraphReader reader;
	Graph airfoil1 = reader.read("input/airfoil1.gi");
#else
	METISGraphReader reader;
	Graph airfoil1 = reader.read("input/airfoil1.graph");
#endif
	double exponent = 2.75;
	BalancedLabelPropagation blp(exponent);
	count k = 4;
	//count n = airfoil1.numberOfNodes();

	// *** airfoil1 graph
	// ML-BLP
	ClusteringGenerator gen;
	Partition partition = gen.makeContinuousBalancedClustering(airfoil1, k);
	std::vector<count> clusterSizes = partition.subsetSizes();
	for (index p = 0; p < k; ++p) {
		DEBUG("size of cluster " , p , ": " , clusterSizes[p]);
	}

	count numVcycles = 50;
	partition = blp.multilevelRun(airfoil1, k);
	for (index vcycle = 1; vcycle < numVcycles; ++vcycle) {
		partition = blp.multilevelRerun(airfoil1, k, partition);
		INFO("cut/balance in cycle " , vcycle , ": " , edgeCut.getQuality(partition, airfoil1) , ", " , GraphClusteringTools::getImbalance(partition));
	}

	double cut = edgeCut.getQuality(partition, airfoil1);
	INFO("BLP number of airfoil1 clusters: " , partition.numberOfSubsets());
	INFO("BLP modularity airfoil1 graph:   " , modularity.getQuality(partition, airfoil1) , ", cut: " , cut);
	clusterSizes = partition.subsetSizes();
	for (index p = 0; p < k; ++p) {
		DEBUG("size of cluster " , p , ": " , clusterSizes[p]);
	}

#if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
	// output to EPS file
	PostscriptWriter psWriter(false);
	psWriter.write(airfoil1, partition, "output/airfoil1-mlblp-4p.eps");
	ClusteringWriter partWriter;
	partWriter.write(partition, "output/airfoil1-mlblp.4p");
#endif
}


// TEST_F(PartitioningAlgoGTest, tryMultilevelKernighanLinOnRealGraph) {
// 	Modularity modularity;
// 	EdgeCut edgeCut;
// #if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
// 	DibapGraphReader reader;
// 	Graph airfoil1 = reader.read("input/airfoil1.gi");
// #else
// 	METISGraphReader reader;
// 	Graph airfoil1 = reader.read("input/airfoil1.graph");
// #endif
// 	KernighanLin partitioner;
// 	count k = 4;

// 	// *** airfoil1 graph
// 	// ML-KL
// 	count numVcycles = 5;
// 	Partition partition = partitioner.multilevelRun(airfoil1, k);
// 	for (index vcycle = 1; vcycle < numVcycles; ++vcycle) {
// 		partition = partitioner.multilevelRerun(airfoil1, k, partition);
// 	}

// 	double cut = edgeCut.getQuality(partition, airfoil1);
// 	INFO("ML-KL number of airfoil1 blocks: " , partition.numberOfSubsets());
// 	INFO("ML-KL modularity airfoil1 graph:   " , modularity.getQuality(partition, airfoil1) , ", cut: " , cut);
// 	std::vector<count> clusterSizes = partition.subsetSizes();
// 	for (index p = 0; p < k; ++p) {
// 		DEBUG("size of cluster " , p , ": " , clusterSizes[p]);
// 	}

// #if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
// 	// output to EPS file
// 	PostscriptWriter psWriter(airfoil1, false);
// 	psWriter.write(partition, "output/airfoil1-mlblp-4p.eps");
// 	ClusteringWriter partWriter;
// 	partWriter.write(partition, "output/airfoil1-mlblp.4p");
// #endif
// }


// TEST_F(PartitioningAlgoGTest, tryKernighanLinStumppOnRealGraph) {
// 	Modularity modularity;
// 	EdgeCut edgeCut;
// #if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
// 	DibapGraphReader reader;
// 	Graph airfoil1 = reader.read("input/airfoil1.gi");
// #else
// 	METISGraphReader reader;
// 	Graph airfoil1 = reader.read("input/airfoil1.graph");
// #endif

// 	KERNIGHAN_LIN partitioner;
// 	count k = 2;

// 	// *** airfoil1 graph
// 	// KL
// 	Partition partition = partitioner.run(airfoil1);

// 	double cut = edgeCut.getQuality(partition, airfoil1);
// 	INFO("ML-KL number of airfoil1 blocks: " , partition.numberOfSubsets());
// 	INFO("ML-KL modularity airfoil1 graph:   " , modularity.getQuality(partition, airfoil1) , ", cut: " , cut);
// 	std::vector<count> clusterSizes = partition.subsetSizes();
// 	for (index p = 0; p < k; ++p) {
// 		DEBUG("size of cluster " , p , ": " , clusterSizes[p]);
// 	}

// #if !defined _WIN32 && !defined _WIN64 && !defined WIN32 && !defined WIN64
// 	// output to EPS file
// 	PostscriptWriter psWriter(airfoil1, false);
// 	psWriter.write(partition, "output/airfoil1-mlblp-4p.eps");
// 	ClusteringWriter partWriter;
// 	partWriter.write(partition, "output/airfoil1-mlblp.4p");
// #endif
// }



} /* namespace NetworKit */

#endif /*NOGTEST */
