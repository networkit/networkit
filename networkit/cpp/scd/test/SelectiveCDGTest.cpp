#include <gtest/gtest.h>
#include <memory>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/community/Conductance.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/SNAPGraphReader.hpp>
#include <networkit/scd/CliqueDetect.hpp>
#include <networkit/scd/GCE.hpp>
#include <networkit/scd/LFMLocal.hpp>
#include <networkit/scd/PageRankNibble.hpp>
#include <networkit/scd/ApproximatePageRank.hpp>
#include <networkit/scd/LocalT.hpp>
#include <networkit/scd/LocalTightnessExpansion.hpp>
#include <networkit/scd/SelectiveCommunityDetector.hpp>
#include <networkit/scd/TCE.hpp>
#include <networkit/scd/TwoPhaseL.hpp>

namespace NetworKit {

class SelectiveCDGTest: public testing::Test {};

TEST_F(SelectiveCDGTest, testRunApproximatePageRank) {
    SNAPGraphReader reader;
    auto G = reader.read("./input/wiki-Vote.txt");

    ApproximatePageRank apr(G, 0.4);
    const auto prVector = apr.run(0);
}

TEST_F(SelectiveCDGTest, testSCD) {
    Aux::Random::setSeed(32, false);
    METISGraphReader reader;
    Graph G = reader.read("input/hep-th.graph");
    // parameters
    node seed = 50;
    std::set<node> seeds = {seed};
    double alpha = 0.1; // loop (or teleport) probability, changed due to DGleich from: // phi * phi / (225.0 * log(100.0 * sqrt(m)));
    double epsilon = 1e-5; // changed due to DGleich from: pow(2, exponent) / (48.0 * B);

    std::vector<std::pair<std::string, std::unique_ptr<SelectiveCommunityDetector>>> algorithms;
    algorithms.emplace_back(std::make_pair(std::string("PageRankNibble"), std::unique_ptr<SelectiveCommunityDetector>(new PageRankNibble(G, alpha, epsilon))));
    algorithms.emplace_back(std::make_pair(std::string("GCE L"), std::unique_ptr<SelectiveCommunityDetector>(new GCE(G, "L"))));
    algorithms.emplace_back(std::make_pair(std::string("GCE M"), std::unique_ptr<SelectiveCommunityDetector>(new GCE(G, "M"))));
    algorithms.emplace_back(std::make_pair(std::string("LFM"), std::unique_ptr<SelectiveCommunityDetector>(new LFMLocal(G, 0.8))));
    algorithms.emplace_back(std::make_pair(std::string("TwoPhaseL"), std::unique_ptr<SelectiveCommunityDetector>(new TwoPhaseL(G))));
    algorithms.emplace_back(std::make_pair(std::string("TCE"), std::unique_ptr<SelectiveCommunityDetector>(new TCE(G))));
    algorithms.emplace_back(std::make_pair(std::string("LTE"), std::unique_ptr<SelectiveCommunityDetector>(new LocalTightnessExpansion(G))));
    algorithms.emplace_back(std::make_pair(std::string("LocalT"), std::unique_ptr<SelectiveCommunityDetector>(new LocalT(G))));


    count idBound = G.upperNodeIdBound();

    for (auto &algIt : algorithms) {
        // run SCD algorithm and partition the graph accordingly
        DEBUG("Call ", algIt.first, "(", seed, ")");
        auto result = algIt.second->run(seeds);
        auto cluster = result[seed];

        // prepare result
        EXPECT_GT(cluster.size(), 0u);
        Partition partition(idBound);
        partition.allToOnePartition();
        partition.toSingleton(seed);
        index id = partition[seed];
        for (auto entry: cluster) {
            partition.moveToSubset(id, entry);
        }

        // evaluate result
        Conductance conductance;
        double targetCond = 1.0;
        double cond = conductance.getQuality(partition, G);
        EXPECT_LT(cond, targetCond);
        INFO("Conductance of ", algIt.first, ": ", cond, "; cluster size: ", cluster.size());
    }
}

TEST_F(SelectiveCDGTest, testGCE) {
    METISGraphReader reader;
    Graph G = reader.read("input/hep-th.graph");

    node seed = 50;
    // Use the "M" score to ensure that the conductance we measure below can only improve
    GCE gce(G, "M");
    auto cluster = gce.expandOneCommunity(seed);

    EXPECT_GT(cluster.size(), 0u);
    double cond1;

    {
        Partition partition(G.upperNodeIdBound());
        partition.allToOnePartition();
        partition.toSingleton(50);
        index id = partition[seed];
        for (auto entry: cluster) {
            partition.moveToSubset(id, entry);
        }

        // evaluate result
        Conductance conductance;
        double targetCond = 1.0;
        cond1 = conductance.getQuality(partition, G);
        EXPECT_LT(cond1, targetCond);
        INFO("Conductance of GCE: ", cond1, "; cluster size: ", cluster.size());
    }

    auto cluster2 = gce.expandOneCommunity(cluster);

    {
        Partition partition(G.upperNodeIdBound());
        partition.allToOnePartition();
        partition.toSingleton(50);
        index id = partition[seed];
        for (auto entry: cluster2) {
            partition.moveToSubset(id, entry);
        }

        // evaluate result
        Conductance conductance;
        // The quality should only improve
        double cond = conductance.getQuality(partition, G);
        EXPECT_LE(cond, cond1);
        INFO("Conductance of GCE2: ", cond, "; cluster size: ", cluster2.size());
    }

    // The cluster should only grow
    for (node u : cluster) {
        EXPECT_NE(cluster2.find(u), cluster2.end());
    }
}

TEST_F(SelectiveCDGTest, testSCDWeighted) {
    Aux::Random::setSeed(23, false);
    METISGraphReader reader;
    Graph G = reader.read("input/lesmis.graph");
    // parameters
    node seed = 20;
    std::set<node> seeds;
    G.forNodes([&](node u) { seeds.insert(u); });
    double alpha = 0.1; // loop (or teleport) probability, changed due to DGleich from: // phi * phi / (225.0 * log(100.0 * sqrt(m)));
    double epsilon = 1e-5; // changed due to DGleich from: pow(2, exponent) / (48.0 * B);

    std::vector<std::pair<std::string, std::unique_ptr<SelectiveCommunityDetector>>> algorithms;
    algorithms.emplace_back(std::make_pair(std::string("PageRankNibble"), std::unique_ptr<SelectiveCommunityDetector>(new PageRankNibble(G, alpha, epsilon))));
    algorithms.emplace_back(std::make_pair(std::string("GCE L"), std::unique_ptr<SelectiveCommunityDetector>(new GCE(G, "L"))));
    algorithms.emplace_back(std::make_pair(std::string("GCE M"), std::unique_ptr<SelectiveCommunityDetector>(new GCE(G, "M"))));
    algorithms.emplace_back(std::make_pair(std::string("LTE"), std::unique_ptr<SelectiveCommunityDetector>(new LocalTightnessExpansion(G))));
    algorithms.emplace_back(std::make_pair(std::string("TCE"), std::unique_ptr<SelectiveCommunityDetector>(new TCE(G))));

    count idBound = G.upperNodeIdBound();

    for (auto &algIt : algorithms) {
        // run SCD algorithm and partition the graph accordingly
        DEBUG("Call ", algIt.first, "(", seed, ")");
        auto result = algIt.second->run(seeds);
        auto cluster = result[seed];

        // prepare result
        EXPECT_GT(cluster.size(), 0u);
        Partition partition(idBound);
        partition.allToOnePartition();
        partition.toSingleton(seed);
        index id = partition[seed];
        for (auto entry: cluster) {
            partition.moveToSubset(id, entry);
        }

        // evaluate result
        Conductance conductance;
        double targetCond = 0.5;
        double cond = conductance.getQuality(partition, G);
        EXPECT_LT(cond, targetCond);
        INFO("Conductance of ", algIt.first, ": ", cond, "; cluster size: ", cluster.size());
    }
}

TEST_F(SelectiveCDGTest, testWeightedCliqueDetect) {
    Graph G(4, true, false);

    G.addEdge(0, 1, 2);
    G.addEdge(0, 2, 1);

    CliqueDetect cliqueDetect(G);

    {
        auto result = cliqueDetect.expandOneCommunity(0);
        EXPECT_EQ(result.size(), 2);
        EXPECT_EQ(result.count(0), 1);
        EXPECT_EQ(result.count(1), 1);
    }

    G.setWeight(0, 2, 3);

    {
        auto result = cliqueDetect.expandOneCommunity(0);
        EXPECT_EQ(result.size(), 2);
        EXPECT_EQ(result.count(0), 1);
        EXPECT_EQ(result.count(2), 1);
    }

    G.setWeight(0, 2, 1);
    G.addEdge(0, 3, 0.5);
    G.addEdge(2, 3, 0.4);

    {
        auto result = cliqueDetect.expandOneCommunity(0);
        EXPECT_EQ(result.size(), 2);
        EXPECT_EQ(result.count(0), 1);
        EXPECT_EQ(result.count(1), 1);
    }

    G.setWeight(2, 3, 0.6);

    {
        auto result = cliqueDetect.expandOneCommunity(0);
        EXPECT_EQ(result.size(), 3);
        EXPECT_EQ(result.count(0), 1);
        EXPECT_EQ(result.count(2), 1);
        EXPECT_EQ(result.count(3), 1);
    }
}

} /* namespace NetworKit */
