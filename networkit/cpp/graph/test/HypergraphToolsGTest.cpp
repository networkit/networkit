/*
 * HypergraphToolsGTest.cpp
 *
 *  Created on: 12.06.2024
 *      Author: Fabian Brandt-Tumescheit
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/graph/HypergraphTools.hpp>

namespace NetworKit {

class HypergraphToolsGTest : public testing::TestWithParam<bool> {
protected:
    bool weighted() const noexcept;
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, HypergraphToolsGTest, testing::Values(false, true));

bool HypergraphToolsGTest::weighted() const noexcept {
    return GetParam();
}

TEST_P(HypergraphToolsGTest, testRandomEdge) {
    constexpr count numNodes = 20;
    constexpr count numEdges = 20;
    constexpr count samples = 100000;
    constexpr double maxAbsoluteError = 0.005;

    auto hGraph = Hypergraph(numNodes, numEdges, weighted());

    std::vector<count> drawCounts(numEdges, 0);
    for (count i = 0; i < samples; i++) {
        ++drawCounts[HypergraphTools::randomEdge(hGraph)];
    }

    for (edgeid eId = 0; eId < numEdges; ++eId) {
        const auto p = static_cast<double>(drawCounts[eId]) / static_cast<double>(samples);
        ASSERT_NEAR(1.0 / numEdges, p, maxAbsoluteError);
    }
}

TEST_P(HypergraphToolsGTest, testRandomEdges) {
    const count numEdges = 20;

    Aux::Random::setSeed(42, false);

    auto hGraph = Hypergraph(0, numEdges, weighted());

    count sampleAll = hGraph.numberOfEdges();
    std::vector<node> allSample = HypergraphTools::randomEdges(hGraph, sampleAll);
    EXPECT_EQ(allSample.size(), sampleAll);
    EXPECT_EQ((std::unordered_set<node>{allSample.begin(), allSample.end()}).size(), sampleAll);

    count sampleMost = 15;
    std::vector<node> mostSample = HypergraphTools::randomEdges(hGraph, sampleMost);
    EXPECT_EQ(mostSample.size(), sampleMost);
    EXPECT_EQ((std::unordered_set<node>{mostSample.begin(), mostSample.end()}).size(), sampleMost);

    count sampleSome = 5;
    std::vector<node> someSample = HypergraphTools::randomEdges(hGraph, sampleSome);
    EXPECT_EQ(someSample.size(), sampleSome);
    EXPECT_EQ((std::unordered_set<node>{someSample.begin(), someSample.end()}).size(), sampleSome);
}

TEST_P(HypergraphToolsGTest, testRandomNode) {
    constexpr count numNodes = 20;
    constexpr count numEdges = 20;
    constexpr count samples = 100000;
    constexpr double maxAbsoluteError = 0.005;

    auto hGraph = Hypergraph(numNodes, numEdges, weighted());

    std::vector<count> drawCounts(numNodes, 0);
    for (count i = 0; i < samples; i++) {
        ++drawCounts[HypergraphTools::randomNode(hGraph)];
    }

    for (node nId = 0; nId < numNodes; ++nId) {
        const auto p = static_cast<double>(drawCounts[nId]) / static_cast<double>(samples);
        ASSERT_NEAR(1.0 / numNodes, p, maxAbsoluteError);
    }
}

TEST_P(HypergraphToolsGTest, testRandomNodes) {
    const count numNodes = 20;

    Aux::Random::setSeed(42, false);

    auto hGraph = Hypergraph(numNodes, 0, weighted());

    count sampleAll = hGraph.numberOfNodes();
    std::vector<node> allSample = HypergraphTools::randomNodes(hGraph, sampleAll);
    EXPECT_EQ(allSample.size(), sampleAll);
    EXPECT_EQ((std::unordered_set<node>{allSample.begin(), allSample.end()}).size(), sampleAll);

    count sampleMost = 15;
    std::vector<node> mostSample = HypergraphTools::randomNodes(hGraph, sampleMost);
    EXPECT_EQ(mostSample.size(), sampleMost);
    EXPECT_EQ((std::unordered_set<node>{mostSample.begin(), mostSample.end()}).size(), sampleMost);

    count sampleSome = 5;
    std::vector<node> someSample = HypergraphTools::randomNodes(hGraph, sampleSome);
    EXPECT_EQ(someSample.size(), sampleSome);
    EXPECT_EQ((std::unordered_set<node>{someSample.begin(), someSample.end()}).size(), sampleSome);
}

TEST_P(HypergraphToolsGTest, testGetIntersection) {
    Hypergraph hGraph(4, 0, weighted());
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({2, 3});

    std::unordered_set<node> intersection0_1{1};
    std::unordered_set<node> intersection1_2{2, 3};
    std::unordered_set<node> intersection2_0{};

    EXPECT_EQ(HypergraphTools::getIntersection(hGraph, 0, 1), intersection0_1);
    EXPECT_EQ(HypergraphTools::getIntersection(hGraph, 1, 2), intersection1_2);
    EXPECT_EQ(HypergraphTools::getIntersection(hGraph, 2, 0), intersection2_0);
}

TEST_P(HypergraphToolsGTest, testGetIntersectionSize) {
    Hypergraph hGraph(4, 0, weighted());
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({2, 3});

    EXPECT_EQ(HypergraphTools::getIntersectionSize(hGraph, 0, 1), 1);
    EXPECT_EQ(HypergraphTools::getIntersectionSize(hGraph, 1, 2), 2);
    EXPECT_EQ(HypergraphTools::getIntersectionSize(hGraph, 2, 0), 0);
}

TEST_P(HypergraphToolsGTest, testCliqueExpansion) {
    Hypergraph hGraph(4, 0, weighted());
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({0, 3});

    Graph cGraph = HypergraphTools::cliqueExpansion(hGraph);

    EXPECT_EQ(cGraph.numberOfNodes(), 4);
    EXPECT_EQ(cGraph.numberOfEdges(), 5);
    ASSERT_TRUE(cGraph.hasEdge(0, 1));
    ASSERT_TRUE(cGraph.hasEdge(1, 2));
    ASSERT_TRUE(cGraph.hasEdge(1, 3));
    ASSERT_TRUE(cGraph.hasEdge(2, 3));
    ASSERT_TRUE(cGraph.hasEdge(0, 3));
}

TEST_P(HypergraphToolsGTest, testLineExpansion) {
    Hypergraph hGraph(4, 0, weighted());
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({0, 3});

    Graph lGraph = HypergraphTools::lineExpansion(hGraph);

    EXPECT_EQ(lGraph.numberOfNodes(), 7);
    EXPECT_EQ(lGraph.numberOfEdges(), 8);

    auto attrNodeRef = lGraph.nodeAttributes().get<node>("node");
    auto attrEdgeRef = lGraph.nodeAttributes().get<edgeid>("edgeid");

    lGraph.forNodes([&](node u) { ASSERT_TRUE(hGraph.hasNode(attrNodeRef[u], attrEdgeRef[u])); });
}

TEST_P(HypergraphToolsGTest, testLineGraph) {
    Hypergraph hGraph(4, 0, weighted());
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({2, 3});

    Graph unweightedLineGraph = HypergraphTools::lineGraph(hGraph, false);
    EXPECT_TRUE(unweightedLineGraph.hasEdge(0, 1));
    EXPECT_TRUE(unweightedLineGraph.hasEdge(1, 2));
    EXPECT_FALSE(unweightedLineGraph.hasEdge(2, 0));
}

TEST_P(HypergraphToolsGTest, testWeightedLineGraph) {
    Hypergraph hGraph(4, 0, weighted());
    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({2, 3});

    Graph weightedLineGraph = HypergraphTools::lineGraph(hGraph, true);
    EXPECT_TRUE(weightedLineGraph.hasEdge(1, 2));
    EXPECT_DOUBLE_EQ(weightedLineGraph.weight(1, 2), 0.5);
    EXPECT_FALSE(weightedLineGraph.hasEdge(2, 0));
    EXPECT_DOUBLE_EQ(weightedLineGraph.weight(2, 0), 0.0);
}

} // namespace NetworKit
