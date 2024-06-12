/*
 * HypergraphGTest.cpp
 *
 *  Created on: 10.06.2024
 *      Author: Fabian Brandt-Tumescheit
 *
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/SimpleHypergraphGenerator.hpp>
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/graph/HypergraphTools.hpp>

namespace NetworKit {

class HypergraphGTest : public testing::TestWithParam<std::tuple<bool>> {
protected:
    bool isHypergraph() const { return !isWeighted(); }
    bool isWeightedHypergraph() const { return isWeighted(); }

    bool isWeighted() const;
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, HypergraphGTest,
                         testing::Values(std::make_tuple(false), std::make_tuple(true)));

bool HypergraphGTest::isWeighted() const {
    return std::get<0>(GetParam());
}

/** NODE MODIFIERS **/

TEST_P(HypergraphGTest, testAddNode) {
    Hypergraph hGraph = Hypergraph();

    ASSERT_FALSE(hGraph.hasNode(0));
    ASSERT_FALSE(hGraph.hasNode(1));
    ASSERT_EQ(0u, hGraph.numberOfNodes());

    hGraph.addNode();
    ASSERT_TRUE(hGraph.hasNode(0));
    ASSERT_FALSE(hGraph.hasNode(1));
    ASSERT_EQ(1u, hGraph.numberOfNodes());

    Hypergraph hGraph2 = Hypergraph(2);
    ASSERT_TRUE(hGraph2.hasNode(0));
    ASSERT_TRUE(hGraph2.hasNode(1));
    ASSERT_FALSE(hGraph2.hasNode(2));
    ASSERT_EQ(2u, hGraph2.numberOfNodes());

    hGraph2.addNode();
    hGraph2.addNode();
    ASSERT_TRUE(hGraph2.hasNode(2));
    ASSERT_TRUE(hGraph2.hasNode(3));
    ASSERT_FALSE(hGraph2.hasNode(4));
    ASSERT_EQ(4u, hGraph2.numberOfNodes());
}

TEST_P(HypergraphGTest, testAddNodes) {
    Hypergraph hGraph = Hypergraph(5);
    hGraph.addNodes(5);

    ASSERT_EQ(hGraph.numberOfNodes(), 10);
    ASSERT_TRUE(hGraph.hasNode(9));

    hGraph.addNodes(90);

    ASSERT_EQ(hGraph.numberOfNodes(), 100);
    ASSERT_TRUE(hGraph.hasNode(99));
}

TEST_P(HypergraphGTest, testRemoveNode) {
    auto testGraph = [&](Hypergraph &hGraph) {
        count numNodes = hGraph.numberOfNodes();
        count maxNodeId = numNodes;
        for (node u = 0; u < maxNodeId; ++u) {
            hGraph.removeNode(u);
            --numNodes;
            ASSERT_EQ(hGraph.numberOfNodes(), numNodes);
            hGraph.forNodes([&](node v) { ASSERT_EQ(hGraph.hasNode(v), v > u); });
        }
    };

    Hypergraph hGraph = SimpleHypergraphGenerator(100, 10, 5).generate();
    testGraph(hGraph);
}

TEST_P(HypergraphGTest, testHasNode) {
    Hypergraph hGraph = Hypergraph(5);

    ASSERT_TRUE(hGraph.hasNode(0));
    ASSERT_TRUE(hGraph.hasNode(1));
    ASSERT_TRUE(hGraph.hasNode(2));
    ASSERT_TRUE(hGraph.hasNode(3));
    ASSERT_TRUE(hGraph.hasNode(4));
    ASSERT_FALSE(hGraph.hasNode(5));
    ASSERT_FALSE(hGraph.hasNode(6));

    hGraph.removeNode(0);
    hGraph.removeNode(2);
    hGraph.addNode();

    ASSERT_FALSE(hGraph.hasNode(0));
    ASSERT_TRUE(hGraph.hasNode(1));
    ASSERT_FALSE(hGraph.hasNode(2));
    ASSERT_TRUE(hGraph.hasNode(3));
    ASSERT_TRUE(hGraph.hasNode(4));
    ASSERT_TRUE(hGraph.hasNode(5));
    ASSERT_FALSE(hGraph.hasNode(6));
}

TEST_P(HypergraphGTest, testRestoreNode) {
    Hypergraph hGraph = Hypergraph(4);

    ASSERT_EQ(4u, hGraph.numberOfNodes());
    ASSERT_TRUE(hGraph.hasNode(0));
    ASSERT_TRUE(hGraph.hasNode(1));
    ASSERT_TRUE(hGraph.hasNode(2));
    ASSERT_TRUE(hGraph.hasNode(3));

    hGraph.removeNode(0);

    ASSERT_EQ(3u, hGraph.numberOfNodes());
    ASSERT_FALSE(hGraph.hasNode(0));
    ASSERT_TRUE(hGraph.hasNode(1));
    ASSERT_TRUE(hGraph.hasNode(2));
    ASSERT_TRUE(hGraph.hasNode(3));

    hGraph.restoreNode(0);

    ASSERT_EQ(4u, hGraph.numberOfNodes());
    ASSERT_TRUE(hGraph.hasNode(0));
    ASSERT_TRUE(hGraph.hasNode(1));
    ASSERT_TRUE(hGraph.hasNode(2));
    ASSERT_TRUE(hGraph.hasNode(3));
}

/** NODE PROPERTIES **/

TEST_P(HypergraphGTest, testDegree) {
    Hypergraph hGraph = SimpleHypergraphGenerator(100, 100, none, false, 10).generate();
    hGraph.forNodes([&](node u) { ASSERT_EQ(10u, hGraph.getDegree(u)); });
}
} // namespace NetworKit
