/*
 * CoarseningGTest.cpp
 *
 *  Created on: 20.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

// networkit-format
#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/coarsening/ClusteringProjector.hpp>
#include <networkit/coarsening/GraphCoarsening.hpp>
#include <networkit/coarsening/MatchingCoarsening.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/community/GraphClusteringTools.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/matching/LocalMaxMatcher.hpp>

namespace NetworKit {

class CoarseningGTest : public testing::Test {};

TEST_F(CoarseningGTest, testDirectedPartitionCoarseningWeights) {
    Graph G(9, true, true);
    for (node u : {1, 3, 4, 6})
        G.removeNode(u);
    G.addEdge(0, 0, 1.0);
    G.addEdge(0, 2, 2.0);
    G.addEdge(2, 0, 3.0);
    G.addEdge(2, 2, 0.5);
    G.addEdge(0, 5, 4.0);
    G.addEdge(2, 5, 1.0);
    G.addEdge(5, 0, 7.0);
    G.addEdge(5, 7, 0.5);
    G.addEdge(7, 5, 1.5);
    G.addEdge(7, 7, 2.0);
    G.addEdge(0, 8, 0.0);
    G.addEdge(2, 8, 0.0);

    Partition partition(G.upperNodeIdBound(), none);
    partition.setUpperBound(18);
    partition[0] = partition[2] = 3;
    partition[5] = partition[7] = 11;
    partition[8] = 17;

    for (bool parallel : {false, true}) {
        ParallelPartitionCoarsening coarsening(G, partition, parallel);
        coarsening.run();
        const Graph &coarse = coarsening.getCoarseGraph();
        const auto &mapping = coarsening.getFineToCoarseNodeMapping();
        ASSERT_EQ(mapping.size(), G.upperNodeIdBound());
        for (node u : {1, 3, 4, 6})
            EXPECT_EQ(mapping[u], none);
        const node a = mapping[0], b = mapping[5], c = mapping[8];
        EXPECT_EQ(mapping[2], a);
        EXPECT_EQ(mapping[7], b);
        EXPECT_NE(a, b);
        EXPECT_NE(a, c);
        EXPECT_NE(b, c);
        ASSERT_TRUE(coarse.isDirected());
        ASSERT_TRUE(coarse.isWeighted());
        ASSERT_TRUE(coarse.checkConsistency());
        EXPECT_EQ(coarse.numberOfNodes(), 3);
        EXPECT_EQ(coarse.numberOfEdges(), 5);
        EXPECT_EQ(coarse.numberOfSelfLoops(), 2);
        EXPECT_DOUBLE_EQ(coarse.weight(a, a), 6.5);
        EXPECT_DOUBLE_EQ(coarse.weight(a, b), 5.0);
        EXPECT_DOUBLE_EQ(coarse.weight(b, a), 7.0);
        EXPECT_DOUBLE_EQ(coarse.weight(b, b), 4.0);
        EXPECT_TRUE(coarse.hasEdge(a, c));
        EXPECT_DOUBLE_EQ(coarse.weight(a, c), 0.0);
        EXPECT_EQ(coarse.degreeIn(c), 1);
        EXPECT_EQ(coarse.degreeOut(c), 0);
        EXPECT_DOUBLE_EQ(coarse.totalEdgeWeight(), 22.5);
        EXPECT_DOUBLE_EQ(coarse.weightedDegree(a, false), 11.5);
        EXPECT_DOUBLE_EQ(coarse.weightedDegreeIn(a), 13.5);
        EXPECT_DOUBLE_EQ(coarse.weightedDegree(b, false), 11.0);
        EXPECT_DOUBLE_EQ(coarse.weightedDegreeIn(b), 9.0);
    }
}

TEST_F(CoarseningGTest, testDirectedPartitionCoarseningObjective) {
    // Independent rational fixture: two asymmetric flow clusters, loops,
    // and a one-way bridge. Q({0,1,2},{3,4,5}) is exactly 36/77.
    Graph G(6, true, true);
    G.addEdge(0, 1, 3.0);
    G.addEdge(1, 0, 1.0);
    G.addEdge(1, 2, 2.0);
    G.addEdge(2, 0, 2.0);
    G.addEdge(2, 3, 0.5);
    G.addEdge(3, 4, 4.0);
    G.addEdge(4, 3, 2.0);
    G.addEdge(4, 5, 1.0);
    G.addEdge(5, 3, 3.0);
    G.addEdge(0, 0, 0.25);
    G.addEdge(5, 5, 0.5);
    Partition partition(std::vector<index>{0, 0, 0, 1, 1, 1});
    Modularity quality;
    ASSERT_NEAR(quality.getQuality(partition, G), 36.0 / 77.0, 1e-14);

    for (bool parallel : {false, true}) {
        ParallelPartitionCoarsening coarsening(G, partition, parallel);
        coarsening.run();
        const Graph &coarse = coarsening.getCoarseGraph();
        const auto &mapping = coarsening.getFineToCoarseNodeMapping();
        ASSERT_TRUE(coarse.isDirected());
        ASSERT_TRUE(coarse.checkConsistency());
        ASSERT_EQ(coarse.numberOfNodes(), 2);
        EXPECT_EQ(coarse.numberOfEdges(), 3);
        EXPECT_DOUBLE_EQ(coarse.totalEdgeWeight(), 77.0 / 4.0);
        EXPECT_DOUBLE_EQ(coarse.weight(mapping[0], mapping[0]), 33.0 / 4.0);
        EXPECT_DOUBLE_EQ(coarse.weight(mapping[0], mapping[3]), 0.5);
        EXPECT_FALSE(coarse.hasEdge(mapping[3], mapping[0]));
        EXPECT_DOUBLE_EQ(coarse.weight(mapping[3], mapping[3]), 21.0 / 2.0);

        Partition separate(2);
        separate.allToSingletons();
        EXPECT_NEAR(quality.getQuality(separate, coarse), 36.0 / 77.0, 1e-14);
        Partition merged(2);
        merged.allToOnePartition();
        EXPECT_NEAR(quality.getQuality(merged, coarse), 0.0, 1e-14);
        ParallelPartitionCoarsening again(coarse, merged, parallel);
        again.run();
        const Graph &one = again.getCoarseGraph();
        ASSERT_TRUE(one.checkConsistency());
        EXPECT_EQ(one.numberOfNodes(), 1);
        EXPECT_EQ(one.numberOfEdges(), 1);
        EXPECT_DOUBLE_EQ(one.weight(0, 0), 77.0 / 4.0);
    }
}

TEST_F(CoarseningGTest, testClusteringProjectorWithOneClustering) {
    ErdosRenyiGenerator gen(100, 0.5);
    Graph G0 = gen.generate();

    // get 1-clustering of G0
    ClusteringGenerator clusteringGen;
    Partition zeta0 = clusteringGen.makeOneClustering(G0);

    // contract G0 according to 1-clusterings
    ParallelPartitionCoarsening contract(G0, zeta0);
    contract.run();
    std::vector<std::vector<node>> maps;
    Graph G1 = contract.getCoarseGraph();
    maps.push_back(contract.getFineToCoarseNodeMapping());

    Partition zeta1 = clusteringGen.makeOneClustering(G1);

    ClusteringProjector project;
    Partition zetaBack = project.projectBackToFinest(zeta1, maps, G0);

    EXPECT_TRUE(GraphClusteringTools::equalClusterings(zeta0, zetaBack, G0))
        << "zeta^{1->0} and zeta^{0} should be identical";
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
    std::vector<std::vector<node>> maps;
    Graph G1 = contract.getCoarseGraph();
    maps.push_back(contract.getFineToCoarseNodeMapping());

    Partition zeta1 = clusteringGen.makeSingletonClustering(G1);

    ClusteringProjector project;
    Partition zetaBack = project.projectBackToFinest(zeta1, maps, G0);

    EXPECT_TRUE(GraphClusteringTools::equalClusterings(zeta0, zetaBack, G0))
        << "zeta^{1->0} and zeta^{0} should be identical";
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

    assert(Gcon.checkConsistency());

    EXPECT_EQ(G.numberOfNodes(), Gcon.numberOfNodes())
        << "graph contracted according to singleton clustering should have the same number of "
           "nodes as original";
    EXPECT_EQ(G.numberOfEdges(), Gcon.numberOfEdges())
        << "graph contracted according to singletons clustering should have the same number of "
           "nodes as original";

    DEBUG("coarsening on random partition");
    count k = 2; // number of clusters in random clustering
    Partition random = clusteringGen.makeRandomClustering(G, k);
    ParallelPartitionCoarsening coarsening2(G, random);
    coarsening2.run();
    Graph GconRand = coarsening2.getCoarseGraph();

    EXPECT_EQ(k, GconRand.numberOfNodes())
        << "graph contracted according to random clustering should have the same number of nodes "
           "as there are clusters.";
    EXPECT_EQ(k + 1, GconRand.numberOfEdges())
        << "graph contracted according to random clustering should have k+1 clusters";
}

TEST_F(CoarseningGTest, testSequentialPartitionCoarseningOnErdosRenyi) {
    count n = 100;

    ErdosRenyiGenerator ERGen(n, 0.5);
    Graph G = ERGen.generate();

    ClusteringGenerator clusteringGen;
    Partition singleton = clusteringGen.makeSingletonClustering(G);

    DEBUG("coarsening on singleton partition");
    ParallelPartitionCoarsening coarsening(G, singleton, false); // parallel by default
    coarsening.run();
    Graph Gcon = coarsening.getCoarseGraph();

    assert(Gcon.checkConsistency());

    EXPECT_EQ(G.numberOfNodes(), Gcon.numberOfNodes())
        << "graph contracted according to singleton clustering should have the same number of "
           "nodes as original";
    EXPECT_EQ(G.numberOfEdges(), Gcon.numberOfEdges())
        << "graph contracted according to singletons clustering should have the same number of "
           "nodes as original";

    DEBUG("coarsening on random partition");
    count k = 2; // number of clusters in random clustering
    Partition random = clusteringGen.makeRandomClustering(G, k);
    ParallelPartitionCoarsening coarsening2(G, random, false);
    coarsening2.run();
    Graph GconRand = coarsening2.getCoarseGraph();

    EXPECT_EQ(k, GconRand.numberOfNodes())
        << "graph contracted according to random clustering should have the same number of nodes "
           "as there are clusters.";
    EXPECT_EQ(k + 1, GconRand.numberOfEdges())
        << "graph contracted according to random clustering should have k+1 clusters";
}

TEST_F(CoarseningGTest, testSequentialPartitionCoarseningEdgeWeights) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    G.addEdge(0, 0); // Add a self loop

    ClusteringGenerator clusteringGen;
    count k = 10; // number of clusters in random clustering
    Partition random = clusteringGen.makeRandomClustering(G, k);

    ParallelPartitionCoarsening coarsening(G, random, false);
    coarsening.run();
    Graph coarseGraph = coarsening.getCoarseGraph();

    ASSERT_EQ(G.totalEdgeWeight(), coarseGraph.totalEdgeWeight());
    std::vector<node> fineToCoarse = coarsening.getFineToCoarseNodeMapping();
    for (const auto &cluster : random.getSubsets()) {
        double degreeSum = 0.0;
        for (node u : cluster) {
            degreeSum += G.weightedDegree(u, true);
        }
        node coarseNode = fineToCoarse[*cluster.begin()];
        ASSERT_EQ(degreeSum, coarseGraph.weightedDegree(coarseNode, true));
    }
}

TEST_F(CoarseningGTest, testParallelPartitionCoarseningEdgeWeights) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");
    G.addEdge(0, 0); // Add a self loop

    ClusteringGenerator clusteringGen;
    count k = 10; // number of clusters in random clustering
    Partition random = clusteringGen.makeRandomClustering(G, k);

    ParallelPartitionCoarsening coarsening(G, random, false);
    coarsening.run();
    Graph coarseGraph = coarsening.getCoarseGraph();

    ASSERT_EQ(G.totalEdgeWeight(), coarseGraph.totalEdgeWeight());
    std::vector<node> fineToCoarse = coarsening.getFineToCoarseNodeMapping();
    for (const auto &cluster : random.getSubsets()) {
        double degreeSum = 0.0;
        for (node u : cluster) {
            degreeSum += G.weightedDegree(u, true);
        }
        node coarseNode = fineToCoarse[*cluster.begin()];
        ASSERT_EQ(degreeSum, coarseGraph.weightedDegree(coarseNode, true));
    }
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

    EXPECT_EQ(Gseq.numberOfEdges(), Gpar.numberOfEdges())
        << "sequential and parallel coarsening should produce the same number of edges";

    auto parMapping = parCoarsening.getFineToCoarseNodeMapping();
    auto seqMapping = seqCoarsening.getFineToCoarseNodeMapping();
    Gseq.forNodes([&](node u) {
        EXPECT_EQ(Gseq.degree(u), Gpar.degree(u)) << "node degrees should be equal for node " << u;
        EXPECT_EQ(Gseq.weightedDegree(u), Gpar.weightedDegree(u))
            << "Weighted degrees should be equald for node " << u;
        EXPECT_EQ(parMapping[u], seqMapping[u]) << "mapping is equal";
    });
}

TEST_F(CoarseningGTest, testPartitionCoarseningOnRealGraph) {
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

    EXPECT_EQ(Gseq.numberOfEdges(), Gpar.numberOfEdges())
        << "sequential and parallel coarsening should produce the same number of edges";

    auto parMapping = parCoarsening.getFineToCoarseNodeMapping();
    auto seqMapping = seqCoarsening.getFineToCoarseNodeMapping();
    Gseq.forNodes([&](node u) {
        EXPECT_EQ(Gseq.degree(u), Gpar.degree(u)) << "node degrees should be equal for node " << u;
        EXPECT_EQ(Gseq.weightedDegree(u), Gpar.weightedDegree(u))
            << "Weighted degrees should be equal for node " << u;
        EXPECT_EQ(parMapping[u], seqMapping[u]) << "mapping is equal";
    });
}

TEST_F(CoarseningGTest, testParallelPartitionCoarseningOnRealGraphWithLoops) {
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

    EXPECT_EQ(Gseq.numberOfEdges(), Gpar.numberOfEdges())
        << "sequential and parallel coarsening should produce the same number of edges";

    auto parMapping = parCoarsening.getFineToCoarseNodeMapping();
    auto seqMapping = seqCoarsening.getFineToCoarseNodeMapping();
    Gseq.forNodes([&](node u) {
        EXPECT_EQ(Gseq.degree(u), Gpar.degree(u)) << "node degrees should be equal for node " << u;
        EXPECT_EQ(Gseq.weightedDegree(u), Gpar.weightedDegree(u))
            << "Weighted degrees should be equal for node " << u;
        EXPECT_EQ(parMapping[u], seqMapping[u]) << "mapping is equal";
    });
}

TEST_F(CoarseningGTest, testMatchingContractor) {
    METISGraphReader reader;
    Graph G = reader.read("input/celegans_metabolic.graph");

    LocalMaxMatcher matcher(G);
    matcher.run();
    Matching matching = matcher.getMatching();
    ASSERT_TRUE(matching.isProper(G));

    MatchingCoarsening coarsener(G, matching);
    coarsener.run();
    Graph coarseG = coarsener.getCoarseGraph();
    std::vector<node> fineToCoarse = coarsener.getFineToCoarseNodeMapping();

    EXPECT_GE(coarseG.numberOfNodes(), 0.5 * G.numberOfNodes());
    EXPECT_EQ(G.totalEdgeWeight(), coarseG.totalEdgeWeight());
    if (G.numberOfEdges() > 0) {
        EXPECT_LT(coarseG.numberOfNodes(), G.numberOfNodes());
    }
}

TEST_F(CoarseningGTest, testMatchingContractorWithSelfLoop) {
    Graph G(2, true);
    G.addEdge(0, 1, 1);
    G.addEdge(0, 0, 2);

    LocalMaxMatcher matcher(G);
    matcher.run();
    Matching matching = matcher.getMatching();
    EXPECT_FALSE(matching.areMatched(0, 0));
    EXPECT_TRUE(matching.areMatched(1, 0));
    ASSERT_TRUE(matching.isProper(G));

    MatchingCoarsening coarsener(G, matching);
    coarsener.run();
    Graph coarseG = coarsener.getCoarseGraph();
    EXPECT_EQ(1u, coarseG.numberOfNodes());

    EXPECT_EQ(G.totalEdgeWeight(), coarseG.totalEdgeWeight());
}

TEST_F(CoarseningGTest, testClusteringProjectorWithNClusterings) {
    ErdosRenyiGenerator gen(100, 0.5);
    Graph Gfine = gen.generate();

    for (int n = 1; n < 10; ++n) {
        // get n-clustering of Gfine
        ClusteringGenerator clusteringGen;
        Partition zetafine = clusteringGen.makeContinuousBalancedClustering(Gfine, n);

        // contract Gfine according to n-clusterings
        ParallelPartitionCoarsening contract(Gfine, zetafine);
        contract.run();
        Graph Gcoarse = contract.getCoarseGraph();

        // revert changes by projecting back into single partition
        Partition zetacoarse = clusteringGen.makeOneClustering(Gcoarse);
        ParallelPartitionCoarsening contractBack(Gcoarse, zetacoarse);
        contractBack.run();
        auto res = contractBack.getCoarseToFineNodeMapping();
        ClusteringProjector CP;
        for (auto const &it : res) {
            for (index i = 0; i < it.second.size(); ++i) {
                auto part = CP.projectBack(Gfine, Gcoarse, it.second, zetacoarse);
                EXPECT_EQ(part.numberOfSubsets(), 1);
                EXPECT_EQ(part.numberOfElements(), n);
            }
        }
    }
}
} /* namespace NetworKit */
