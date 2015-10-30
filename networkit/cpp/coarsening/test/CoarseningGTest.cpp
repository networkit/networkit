/*
 * CoarseningGTest.cpp
 *
 *  Created on: 20.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "CoarseningGTest.h"

#include "../../auxiliary/Log.h"
#include "../../community/ClusteringGenerator.h"
#include "../../coarsening/ClusteringProjector.h"
#include "../../community/GraphClusteringTools.h"
#include "../../generators/ErdosRenyiGenerator.h"
#include "../../coarsening/ParallelPartitionCoarsening.h"
#include "../../io/METISGraphReader.h"

namespace NetworKit {

TEST_F(CoarseningGTest, testClusteringProjectorWithOneClustering) {
	ErdosRenyiGenerator gen(100, 0.5);
	Graph G0 = gen.generate();

	// get 1-clustering of G0
	ClusteringGenerator clusteringGen;
	Partition zeta0 = clusteringGen.makeOneClustering(G0);

	// contract G0 according to 1-clusterings
	ParallelPartitionCoarsening contract(G0,zeta0);
	contract.run();
	std::vector<std::vector<node> > maps;
	Graph G1 = contract.getCoarseGraph();
	maps.push_back(contract.getNodeMapping());

	Partition zeta1 = clusteringGen.makeOneClustering(G1);

	ClusteringProjector project;
	Partition zetaBack = project.projectBackToFinest(zeta1, maps, G0);

	EXPECT_TRUE(GraphClusteringTools::equalClusterings(zeta0, zetaBack, G0)) << "zeta^{1->0} and zeta^{0} should be identical";
}


TEST_F(CoarseningGTest, testClusteringProjectorWithSingletonClustering) {
	ErdosRenyiGenerator gen(100, 0.5);
	Graph G0 = gen.generate();

	// get 1-clustering of G0
	ClusteringGenerator clusteringGen;
	Partition zeta0 = clusteringGen.makeSingletonClustering(G0);

	// contract G0 according to 1-clusterings
	ParallelPartitionCoarsening contract(G0, zeta0);
	contract.run();
	std::vector<std::vector<node> > maps;
	Graph G1 = contract.getCoarseGraph();
	maps.push_back(contract.getNodeMapping());

	Partition zeta1 = clusteringGen.makeSingletonClustering(G1);

	ClusteringProjector project;
	Partition zetaBack = project.projectBackToFinest(zeta1, maps, G0);

	EXPECT_TRUE(GraphClusteringTools::equalClusterings(zeta0, zetaBack, G0)) << "zeta^{1->0} and zeta^{0} should be identical";
}


TEST_F(CoarseningGTest, testParallelPartitionCoarseningOnErdosRenyi) {
	count n = 100;

	ErdosRenyiGenerator ERGen(n, 0.5);
	Graph G = ERGen.generate();

	ClusteringGenerator clusteringGen;
	Partition singleton = clusteringGen.makeSingletonClustering(G);


	DEBUG("coarsening on singleton partition");
	ParallelPartitionCoarsening coarsening(G, singleton);
	coarsening.run();
	Graph Gcon = coarsening.getCoarseGraph();

	assert (Gcon.checkConsistency());

	EXPECT_EQ(G.numberOfNodes(), Gcon.numberOfNodes())
			<< "graph contracted according to singleton clustering should have the same number of nodes as original";
	EXPECT_EQ(G.numberOfEdges(), Gcon.numberOfEdges())
			<< "graph contracted according to singletons clustering should have the same number of nodes as original";


	DEBUG("coarsening on random partition");
	count k = 2; // number of clusters in random clustering
	Partition random = clusteringGen.makeRandomClustering(G, k);
	ParallelPartitionCoarsening coarsening2(G, random);
	coarsening2.run();
	Graph GconRand = coarsening2.getCoarseGraph();

	EXPECT_EQ(k, GconRand.numberOfNodes())
			<< "graph contracted according to random clustering should have the same number of nodes as there are clusters.";
	EXPECT_EQ(k + 1, GconRand.numberOfEdges()) << "graph contracted according to random clustering should have k+1 clusters";
}

TEST_F(CoarseningGTest, testParallelPartitionCoarseningOnErdosRenyiWithGraphBuilder) {
	count n = 100;

	ErdosRenyiGenerator ERGen(n, 0.5);
	Graph G = ERGen.generate();

	ClusteringGenerator clusteringGen;
	Partition singleton = clusteringGen.makeSingletonClustering(G);


	DEBUG("coarsening on singleton partition");
	ParallelPartitionCoarsening coarsening(G, singleton); // uses graph builder by default
	coarsening.run();
	Graph Gcon = coarsening.getCoarseGraph();

	assert (Gcon.checkConsistency());

	EXPECT_EQ(G.numberOfNodes(), Gcon.numberOfNodes())
			<< "graph contracted according to singleton clustering should have the same number of nodes as original";
	EXPECT_EQ(G.numberOfEdges(), Gcon.numberOfEdges())
			<< "graph contracted according to singletons clustering should have the same number of nodes as original";

	DEBUG("coarsening on random partition");
	count k = 2; // number of clusters in random clustering
	Partition random = clusteringGen.makeRandomClustering(G, k);
	ParallelPartitionCoarsening coarsening2(G, random);
	coarsening2.run();
	Graph GconRand = coarsening2.getCoarseGraph();

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

	ParallelPartitionCoarsening parCoarsening(G, random);
	parCoarsening.run();

	ParallelPartitionCoarsening seqCoarsening(G, random, false);
	seqCoarsening.run();

	Graph Gpar = parCoarsening.getCoarseGraph();
	EXPECT_EQ(k, Gpar.numberOfNodes());

	Graph Gseq = seqCoarsening.getCoarseGraph();
	EXPECT_EQ(k, Gseq.numberOfNodes());

	EXPECT_EQ(Gseq.numberOfEdges(), Gpar.numberOfEdges()) << "sequential and parallel coarsening should produce the same number of edges";

	auto parMapping = parCoarsening.getNodeMapping();
	auto seqMapping = seqCoarsening.getNodeMapping();
	Gseq.forNodes([&](node u) {
		EXPECT_EQ(Gseq.degree(u), Gpar.degree(u)) << "node degrees should be equal for node " << u;
		EXPECT_EQ(Gseq.weightedDegree(u), Gpar.weightedDegree(u)) << "Weighted degrees should be equald for node " << u;
		EXPECT_EQ(parMapping[u], seqMapping[u]) << "mapping is equal";
	});

}

TEST_F(CoarseningGTest, testParallelPartitionCoarseningOnRealGraphWithGraphBuilder) {
	METISGraphReader reader;
	Graph G = reader.read("input/celegans_metabolic.graph");

	ClusteringGenerator clusteringGen;
	count k = 10; // number of clusters in random clustering
	Partition random = clusteringGen.makeRandomClustering(G, k);

	ParallelPartitionCoarsening parCoarsening(G, random, true);
	parCoarsening.run();
	Graph Gpar = parCoarsening.getCoarseGraph();
	EXPECT_EQ(random.numberOfSubsets(), Gpar.numberOfNodes());

	ParallelPartitionCoarsening seqCoarsening(G, random, false);
	seqCoarsening.run();
	Graph Gseq = seqCoarsening.getCoarseGraph();
	EXPECT_EQ(random.numberOfSubsets(), Gseq.numberOfNodes());

	EXPECT_EQ(Gseq.numberOfEdges(), Gpar.numberOfEdges()) << "sequential and parallel coarsening should produce the same number of edges";

	auto parMapping = parCoarsening.getNodeMapping();
	auto seqMapping = seqCoarsening.getNodeMapping();
	Gseq.forNodes([&](node u) {
		EXPECT_EQ(Gseq.degree(u), Gpar.degree(u)) << "node degrees should be equal for node " << u;
		EXPECT_EQ(Gseq.weightedDegree(u), Gpar.weightedDegree(u)) << "Weighted degrees should be equal for node " << u;
		EXPECT_EQ(parMapping[u], seqMapping[u]) << "mapping is equal";
	});
}

TEST_F(CoarseningGTest, testParallelPartitionCoarseningOnRealGraphWithGraphBuilderAndLoops) {
	METISGraphReader reader;
	Graph G = reader.read("input/celegans_metabolic.graph");
	G.addEdge(0, 0);

	ClusteringGenerator clusteringGen;
	count k = 10; // number of clusters in random clustering
	Partition random = clusteringGen.makeRandomClustering(G, k);

	ParallelPartitionCoarsening parCoarsening(G, random, true);
	parCoarsening.run();
	Graph Gpar = parCoarsening.getCoarseGraph();
	EXPECT_EQ(random.numberOfSubsets(), Gpar.numberOfNodes());

	ParallelPartitionCoarsening seqCoarsening(G, random, false);
	seqCoarsening.run();
	Graph Gseq = seqCoarsening.getCoarseGraph();
	EXPECT_EQ(random.numberOfSubsets(), Gseq.numberOfNodes());

	EXPECT_EQ(Gseq.numberOfEdges(), Gpar.numberOfEdges()) << "sequential and parallel coarsening should produce the same number of edges";

	auto parMapping = parCoarsening.getNodeMapping();
	auto seqMapping = seqCoarsening.getNodeMapping();
	Gseq.forNodes([&](node u) {
		EXPECT_EQ(Gseq.degree(u), Gpar.degree(u)) << "node degrees should be equal for node " << u;
		EXPECT_EQ(Gseq.weightedDegree(u), Gpar.weightedDegree(u)) << "Weighted degrees should be equal for node " << u;
		EXPECT_EQ(parMapping[u], seqMapping[u]) << "mapping is equal";
	});
}

} /* namespace NetworKit */

#endif /*NOGTEST */
