#include <gtest/gtest.h>
#include <memory>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/community/Conductance.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/SNAPGraphReader.hpp>
#include <networkit/scd/GCE.hpp>
#include <networkit/scd/PageRankNibble.hpp>
#include <networkit/scd/ApproximatePageRank.hpp>
#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {

class SCDGTest2: public testing::Test {};

TEST_F(SCDGTest2, testRunApproximatePageRank) {
    SNAPGraphReader reader;
    auto G = reader.read("./input/wiki-Vote.txt");

    ApproximatePageRank apr(G, 0.4);
    const auto prVector= apr.run(0);
}

TEST_F(SCDGTest2, testSCD) {
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
        partition.toSingleton(50);
        index id = partition[seed];
        for (auto entry: cluster) {
            partition.moveToSubset(id, entry);
        }

        // evaluate result
        Conductance conductance;
        double targetCond = 0.4;
        double cond = conductance.getQuality(partition, G);
        EXPECT_LT(cond, targetCond);
        INFO("Conductance of ", algIt.first, ": ", cond, "; cluster size: ", cluster.size());
    }
}


} /* namespace NetworKit */
