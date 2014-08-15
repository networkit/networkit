/*
 * CoarseningGTest.cpp
 *
 *  Created on: 20.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "CoarseningGTest.h"

#include "../../auxiliary/Log.h"
#include "../../graph/GraphGenerator.h"
#include "../../community/ClusteringGenerator.h"
#include "../../coarsening/ClusterContractor.h"
#include "../../coarsening/PartitionCoarsening.h"
#include "../../coarsening/ClusteringProjector.h"
#include "../../community/GraphClusteringTools.h"
#include "../../generators/ErdosRenyiGenerator.h"
#include "../../coarsening/ParallelPartitionCoarsening.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {

TEST_F(CoarseningGTest, testClusterContractor) {
	GraphGenerator graphGen;
	int64_t n = 100;
	Graph G = graphGen.makeErdosRenyiGraph(n, 0.5);

	ClusteringGenerator clusteringGen;
	Partition singleton = clusteringGen.makeSingletonClustering(G);


	ClusterContractor contracter;
	auto conSingletonPair = contracter.run(G, singleton);
	Graph Gcon = conSingletonPair.first;

	EXPECT_EQ(G.numberOfNodes(), Gcon.numberOfNodes())
			<< "graph contracted according to singleton clustering should have the same number of nodes as original";
	EXPECT_EQ(G.numberOfEdges(), Gcon.numberOfEdges())
			<< "graph contracted according to singletons clustering should have the same number of nodes as original";

	count k = 2; // number of clusters in random clustering
	Partition random = clusteringGen.makeRandomClustering(G, k);
	auto conRandPair = contracter.run(G, random);
	Graph GconRand = conRandPair.first;

	EXPECT_EQ(k, GconRand.numberOfNodes())
			<< "graph contracted according to random clustering should have the same number of nodes as there are clusters.";

}

TEST_F(CoarseningGTest, testPartitionCoarsening) {
	GraphGenerator graphGen;
	count n = 100;
	Graph G = graphGen.makeErdosRenyiGraph(n, 0.5);

	ClusteringGenerator clusteringGen;
	Partition singleton = clusteringGen.makeSingletonClustering(G);


	PartitionCoarsening coarsening;
	auto conSingletonPair = coarsening.run(G, singleton);
	Graph Gcon = conSingletonPair.first;

	EXPECT_EQ(G.numberOfNodes(), Gcon.numberOfNodes())
			<< "graph contracted according to singleton clustering should have the same number of nodes as original";
	EXPECT_EQ(G.numberOfEdges(), Gcon.numberOfEdges())
			<< "graph contracted according to singletons clustering should have the same number of nodes as original";

	count k = 2; // number of clusters in random clustering
	Partition random = clusteringGen.makeRandomClustering(G, k);
	auto conRandPair = coarsening.run(G, random);
	Graph GconRand = conRandPair.first;

	EXPECT_EQ(k, GconRand.numberOfNodes())
			<< "graph contracted according to random clustering should have the same number of nodes as there are clusters.";

}

TEST_F(CoarseningGTest, testClusteringProjectorWithOneClustering) {
	GraphGenerator graphGen;
	int64_t n = 100;
	Graph G0 = graphGen.makeErdosRenyiGraph(n, 0.5);

	// get 1-clustering of G0
	ClusteringGenerator clusteringGen;
	Partition zeta0 = clusteringGen.makeOneClustering(G0);

	// contract G0 according to 1-clusterings
	ClusterContractor contract;
	auto con = contract.run(G0, zeta0);
	std::vector<std::vector<node> > maps;
	Graph G1 = con.first;
	maps.push_back(con.second);

	Partition zeta1 = clusteringGen.makeOneClustering(G1);

	ClusteringProjector project;
	Partition zetaBack = project.projectBackToFinest(zeta1, maps, G0);

	EXPECT_TRUE(GraphClusteringTools::equalClusterings(zeta0, zetaBack, G0)) << "zeta^{1->0} and zeta^{0} should be identical";
}


TEST_F(CoarseningGTest, testClusteringProjectorWithSingletonClustering) {
	GraphGenerator graphGen;
	int64_t n = 100;
	Graph G0 = graphGen.makeErdosRenyiGraph(n, 0.5);

	// get 1-clustering of G0
	ClusteringGenerator clusteringGen;
	Partition zeta0 = clusteringGen.makeSingletonClustering(G0);

	// contract G0 according to 1-clusterings
	ClusterContractor contract;
	auto con = contract.run(G0, zeta0);
	std::vector<std::vector<node> > maps;
	Graph G1 = con.first;
	maps.push_back(con.second);

	Partition zeta1 = clusteringGen.makeSingletonClustering(G1);

	ClusteringProjector project;
	Partition zetaBack = project.projectBackToFinest(zeta1, maps, G0);

	EXPECT_TRUE(GraphClusteringTools::equalClusterings(zeta0, zetaBack, G0)) << "zeta^{1->0} and zeta^{0} should be identical";
}


TEST_F(CoarseningGTest, testParallelPartitionCoarsening) {
	count n = 100;

	ErdosRenyiGenerator ERGen(n, 0.5);
	Graph G = ERGen.generate();

	ClusteringGenerator clusteringGen;
	Partition singleton = clusteringGen.makeSingletonClustering(G);


	DEBUG("coarsening on singleton partition");
	ParallelPartitionCoarsening coarsening;
	auto conSingletonPair = coarsening.run(G, singleton);
	Graph Gcon = conSingletonPair.first;

	assert (Gcon.consistencyCheck());

	EXPECT_EQ(G.numberOfNodes(), Gcon.numberOfNodes())
			<< "graph contracted according to singleton clustering should have the same number of nodes as original";
	EXPECT_EQ(G.numberOfEdges(), Gcon.numberOfEdges())
			<< "graph contracted according to singletons clustering should have the same number of nodes as original";


	DEBUG("coarsening on random partition");
	count k = 2; // number of clusters in random clustering
	Partition random = clusteringGen.makeRandomClustering(G, k);
	auto conRandPair = coarsening.run(G, random);
	Graph GconRand = conRandPair.first;

	EXPECT_EQ(k, GconRand.numberOfNodes())
			<< "graph contracted according to random clustering should have the same number of nodes as there are clusters.";
	EXPECT_EQ(k + 1, GconRand.numberOfEdges()) << "graph contracted according to random clustering should have k+1 clusters";

}

TEST_F(CoarseningGTest, testParallelPartitionCoarseningOnRealGraph) {
	METISGraphReader reader;
	Graph G = reader.read("input/celegans_metabolic.graph");


	ClusteringGenerator clusteringGen;
	count k = 10; // number of clusters in random clustering
	Partition random = clusteringGen.makeRandomClustering(G, k);

	ParallelPartitionCoarsening parCoarsening;
	auto parResult = parCoarsening.run(G, random);

	ClusterContractor seqCoarsening;
	auto seqResult = seqCoarsening.run(G, random);

	Graph Gpar = parResult.first;
	EXPECT_EQ(k, Gpar.numberOfNodes());

	Graph Gseq = seqResult.first;
	EXPECT_EQ(k, Gseq.numberOfNodes());

	EXPECT_EQ(Gseq.numberOfEdges(), Gpar.numberOfEdges()) << "sequential and parallel coarsening should produce the same number of edges";

	Gseq.forNodes([&](node u){
		EXPECT_EQ(Gseq.degree(u), Gpar.degree(u)) << "node degrees should be equal";
		EXPECT_EQ(parResult.second[u], seqResult.second[u]) << "mapping is equal";
	});

}


} /* namespace NetworKit */

#endif /*NOGTEST */
