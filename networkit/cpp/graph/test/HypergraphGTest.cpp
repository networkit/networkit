/*
 * HypergraphGTest.cpp
 *
 *  Created on: 10.06.2024
 *      Author: Fabian Brandt-Tumescheit
 *
 */

#include <gmock/gmock.h>
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

/** BASIC PROPERTIES **/

TEST_P(HypergraphGTest, testEmptyHypergraph) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());

    EXPECT_TRUE(hGraph.isEmpty());

    hGraph.addNode();

    EXPECT_FALSE(hGraph.isEmpty());
}

TEST_P(HypergraphGTest, testNumNodes) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());

    EXPECT_EQ(0u, hGraph.numberOfNodes());

    hGraph.addNode();

    EXPECT_EQ(1u, hGraph.numberOfNodes());
}

TEST_P(HypergraphGTest, testNumEdges) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());

    EXPECT_EQ(0u, hGraph.numberOfEdges());

    hGraph.addEdge();

    EXPECT_EQ(1u, hGraph.numberOfEdges());
}

TEST_P(HypergraphGTest, testUpperNodeIdBound) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());

    EXPECT_EQ(0u, hGraph.upperNodeIdBound());

    hGraph.addNodes(10);

    EXPECT_EQ(10u, hGraph.upperNodeIdBound());
}

TEST_P(HypergraphGTest, testUpperEdgeIdBound) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());

    EXPECT_EQ(0u, hGraph.upperEdgeIdBound());

    hGraph.addEdge();
    hGraph.addEdge();
    hGraph.addEdge();

    EXPECT_EQ(3u, hGraph.upperEdgeIdBound());
}

/** NODE MODIFIERS **/

TEST_P(HypergraphGTest, testAddNode) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());

    EXPECT_FALSE(hGraph.hasNode(0));
    EXPECT_FALSE(hGraph.hasNode(1));
    EXPECT_EQ(0u, hGraph.numberOfNodes());

    hGraph.addNode();
    EXPECT_TRUE(hGraph.hasNode(0));
    EXPECT_FALSE(hGraph.hasNode(1));
    EXPECT_EQ(1u, hGraph.numberOfNodes());

    Hypergraph hGraph2 = Hypergraph(2, 0, isWeightedHypergraph());
    EXPECT_TRUE(hGraph2.hasNode(0));
    EXPECT_TRUE(hGraph2.hasNode(1));
    EXPECT_FALSE(hGraph2.hasNode(2));
    EXPECT_EQ(2u, hGraph2.numberOfNodes());

    hGraph2.addNode();
    hGraph2.addNode();
    EXPECT_TRUE(hGraph2.hasNode(2));
    EXPECT_TRUE(hGraph2.hasNode(3));
    EXPECT_FALSE(hGraph2.hasNode(4));
    EXPECT_EQ(4u, hGraph2.numberOfNodes());
}

TEST_P(HypergraphGTest, testAddNodes) {
    Hypergraph hGraph = Hypergraph(5, 0, isWeightedHypergraph());
    hGraph.addNodes(5);

    EXPECT_EQ(hGraph.numberOfNodes(), 10);
    EXPECT_TRUE(hGraph.hasNode(9));

    hGraph.addNodes(90);

    EXPECT_EQ(hGraph.numberOfNodes(), 100);
    EXPECT_TRUE(hGraph.hasNode(99));
}

TEST_P(HypergraphGTest, testAddNodeToNonExisting) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());

    EXPECT_THROW(hGraph.addNodeTo({0}, 0), std::runtime_error);
}

TEST_P(HypergraphGTest, testAddNodeToWithoutNodeId) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());
    hGraph.addEdge();

    hGraph.addNodeTo({0});

    EXPECT_EQ(1u, hGraph.numberOfNodes());
    EXPECT_TRUE(hGraph.hasNode(0));
}

TEST_P(HypergraphGTest, testAddNodeToWithNodeId) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());
    hGraph.addEdge();
    hGraph.addNode();

    hGraph.addNodeTo({0}, 0);

    EXPECT_EQ(1u, hGraph.numberOfNodes());
    EXPECT_TRUE(hGraph.hasNode(0));
}

TEST_P(HypergraphGTest, testAddNodesToNonExisting) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());

    EXPECT_THROW(hGraph.addNodesTo({0}, 0), std::runtime_error);
}

TEST_P(HypergraphGTest, testAddNodesToWithoutEdgeId) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());
    hGraph.addNode();

    hGraph.addNodesTo({0});

    EXPECT_EQ(1u, hGraph.numberOfEdges());
    EXPECT_TRUE(hGraph.hasEdge(0));
}

TEST_P(HypergraphGTest, testAddNodesToWithEdgeId) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());
    hGraph.addEdge();
    hGraph.addNode();

    hGraph.addNodeTo({0}, 0);

    EXPECT_TRUE(hGraph.hasNode(0, 0));
}

TEST_F(HypergraphGTest, testRemoveNode) {
    auto testGraph = [&](Hypergraph &hGraph) {
        count numNodes = hGraph.numberOfNodes();
        count maxNodeId = numNodes;
        for (node u = 0; u < maxNodeId; ++u) {
            hGraph.removeNode(u);
            --numNodes;
            EXPECT_EQ(hGraph.numberOfNodes(), numNodes);
            hGraph.forNodes([&](node v) { EXPECT_EQ(hGraph.hasNode(v), v > u); });
        }
    };

    Hypergraph hGraph = SimpleHypergraphGenerator(100, 10, 5).generate();
    testGraph(hGraph);
}

TEST_P(HypergraphGTest, testRemoveNodeFrom) {
    Hypergraph hGraph = Hypergraph(4, 0, isWeightedHypergraph());

    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({0, 3});

    hGraph.removeNodeFrom(0, 0);
    hGraph.removeNodeFrom(2, 1);
    hGraph.removeNodeFrom(3, 2);

    EXPECT_FALSE(hGraph.hasNode(0, 0));
    EXPECT_TRUE(hGraph.hasNode(1, 0));
    EXPECT_TRUE(hGraph.hasNode(1, 1));
    EXPECT_FALSE(hGraph.hasNode(2, 1));
    EXPECT_TRUE(hGraph.hasNode(3, 1));
    EXPECT_TRUE(hGraph.hasNode(0, 2));
    EXPECT_FALSE(hGraph.hasNode(3, 2));
}

TEST_P(HypergraphGTest, testHasNode) {
    Hypergraph hGraph = Hypergraph(5, 0, isWeightedHypergraph());

    EXPECT_TRUE(hGraph.hasNode(0));
    EXPECT_TRUE(hGraph.hasNode(1));
    EXPECT_TRUE(hGraph.hasNode(2));
    EXPECT_TRUE(hGraph.hasNode(3));
    EXPECT_TRUE(hGraph.hasNode(4));
    EXPECT_FALSE(hGraph.hasNode(5));
    EXPECT_FALSE(hGraph.hasNode(6));

    hGraph.removeNode(0);
    hGraph.removeNode(2);
    hGraph.addNode();

    EXPECT_FALSE(hGraph.hasNode(0));
    EXPECT_TRUE(hGraph.hasNode(1));
    EXPECT_FALSE(hGraph.hasNode(2));
    EXPECT_TRUE(hGraph.hasNode(3));
    EXPECT_TRUE(hGraph.hasNode(4));
    EXPECT_TRUE(hGraph.hasNode(5));
    EXPECT_FALSE(hGraph.hasNode(6));
}

TEST_P(HypergraphGTest, testHasNodeInHyperedge) {
    Hypergraph hGraph = Hypergraph(4, 0, isWeightedHypergraph());

    hGraph.addEdge({0, 1});
    hGraph.addEdge({1, 2, 3});
    hGraph.addEdge({0, 3});

    EXPECT_TRUE(hGraph.hasNode(0, 0));
    EXPECT_TRUE(hGraph.hasNode(1, 0));
    EXPECT_TRUE(hGraph.hasNode(2, 1));
    EXPECT_TRUE(hGraph.hasNode(3, 1));
    EXPECT_TRUE(hGraph.hasNode(3, 2));
    EXPECT_FALSE(hGraph.hasNode(2, 0));
    EXPECT_FALSE(hGraph.hasNode(2, 2));
    EXPECT_FALSE(hGraph.hasNode(3, 0));

    hGraph.removeNode(2);

    EXPECT_FALSE(hGraph.hasNode(2, 1));
}

TEST_P(HypergraphGTest, testRestoreNode) {
    Hypergraph hGraph = Hypergraph(4, 0, isWeightedHypergraph());

    EXPECT_EQ(4u, hGraph.numberOfNodes());
    EXPECT_TRUE(hGraph.hasNode(0));
    EXPECT_TRUE(hGraph.hasNode(1));
    EXPECT_TRUE(hGraph.hasNode(2));
    EXPECT_TRUE(hGraph.hasNode(3));

    hGraph.removeNode(0);

    EXPECT_EQ(3u, hGraph.numberOfNodes());
    EXPECT_FALSE(hGraph.hasNode(0));
    EXPECT_TRUE(hGraph.hasNode(1));
    EXPECT_TRUE(hGraph.hasNode(2));
    EXPECT_TRUE(hGraph.hasNode(3));

    hGraph.restoreNode(0);

    EXPECT_EQ(4u, hGraph.numberOfNodes());
    EXPECT_TRUE(hGraph.hasNode(0));
    EXPECT_TRUE(hGraph.hasNode(1));
    EXPECT_TRUE(hGraph.hasNode(2));
    EXPECT_TRUE(hGraph.hasNode(3));
}

/** NODE PROPERTIES **/

TEST_F(HypergraphGTest, testDegree) {
    Hypergraph hGraph = SimpleHypergraphGenerator(100, 100, none, false, 10).generate();
    hGraph.forNodes([&](node u) { EXPECT_EQ(10u, hGraph.degree(u)); });
}

TEST_P(HypergraphGTest, testWeightedDegree) {
    Hypergraph hGraph = Hypergraph(1, 3, isWeightedHypergraph());
    hGraph.addNodeTo({0}, 0);
    hGraph.addNodeTo({1}, 0);
    hGraph.addNodeTo({2}, 0);

    hGraph.forEdges(
        [&](edgeid eid) { hGraph.setEdgeWeight(eid, static_cast<edgeweight>(eid + 1)); });

    EXPECT_DOUBLE_EQ(hGraph.weightedDegree(0), isWeightedHypergraph() ? 6.0 : 3.0);
}

TEST_P(HypergraphGTest, testNeighborsWithoutDuplicates) {
    Hypergraph hGraph = Hypergraph(3, 3, isWeightedHypergraph());
    hGraph.addNodesTo({0, 1}, 0);
    hGraph.addNodesTo({1}, 1);
    hGraph.addNodesTo({1, 2}, 2);

    auto neighbors = hGraph.getNeighbors(1);

    EXPECT_THAT(neighbors, testing::UnorderedElementsAre(0, 2));
}

TEST_P(HypergraphGTest, testNeighborsWithDuplicates) {
    Hypergraph hGraph = Hypergraph(3, 3, isWeightedHypergraph());
    hGraph.addNodesTo({0, 1}, 0);
    hGraph.addNodesTo({1}, 1);
    hGraph.addNodesTo({0, 1}, 2);

    auto neighbors = hGraph.getNeighbors(1);

    EXPECT_EQ(neighbors.size(), 1);
    EXPECT_THAT(neighbors, testing::ElementsAre(0));
}

TEST_P(HypergraphGTest, testGetSetNodeWeight) {
    Hypergraph hGraph = Hypergraph(3, 0, isWeightedHypergraph());

    EXPECT_DOUBLE_EQ(hGraph.getNodeWeight(0), defaultNodeWeight);

    hGraph.setNodeWeight(1, 2.0);

    EXPECT_DOUBLE_EQ(hGraph.getNodeWeight(1), isWeightedHypergraph() ? 2.0 : defaultNodeWeight);

    hGraph.setNodeWeight(2, -2.0);

    EXPECT_DOUBLE_EQ(hGraph.getNodeWeight(2), isWeightedHypergraph() ? -2.0 : defaultNodeWeight);
}

/** EDGE MODIFIERS **/

TEST_P(HypergraphGTest, testaddEdge) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());

    EXPECT_FALSE(hGraph.hasEdge(0));
    EXPECT_FALSE(hGraph.hasEdge(1));
    EXPECT_EQ(0u, hGraph.numberOfEdges());

    hGraph.addEdge();
    EXPECT_TRUE(hGraph.hasEdge(0));
    EXPECT_FALSE(hGraph.hasEdge(1));
    EXPECT_EQ(1u, hGraph.numberOfEdges());

    Hypergraph hGraph2 = Hypergraph(0, 2, isWeightedHypergraph());
    EXPECT_TRUE(hGraph2.hasEdge(0));
    EXPECT_TRUE(hGraph2.hasEdge(1));
    EXPECT_FALSE(hGraph2.hasEdge(2));
    EXPECT_EQ(2u, hGraph2.numberOfEdges());

    hGraph2.addEdge();
    hGraph2.addEdge();
    EXPECT_TRUE(hGraph2.hasEdge(2));
    EXPECT_TRUE(hGraph2.hasEdge(3));
    EXPECT_FALSE(hGraph2.hasEdge(4));
    EXPECT_EQ(4u, hGraph2.numberOfEdges());
}

TEST_P(HypergraphGTest, testaddEdgeWithNodes) {
    Hypergraph hGraph = Hypergraph(2, 0, isWeightedHypergraph());

    hGraph.addEdge({0, 1});
    EXPECT_TRUE(hGraph.hasEdge(0));
    EXPECT_EQ(1u, hGraph.numberOfEdges());
}

TEST_P(HypergraphGTest, testaddEdgeWithNodesAddMissing) {
    Hypergraph hGraph = Hypergraph(0, 0, isWeightedHypergraph());

    hGraph.addEdge({0, 1}, true);
    EXPECT_TRUE(hGraph.hasEdge(0));
    EXPECT_TRUE(hGraph.hasNode(0, 0));
    EXPECT_TRUE(hGraph.hasNode(1, 0));
    EXPECT_EQ(1u, hGraph.numberOfEdges());
}

TEST_P(HypergraphGTest, testHasEdge) {
    Hypergraph hGraph = Hypergraph(0, 5, isWeightedHypergraph());

    EXPECT_TRUE(hGraph.hasEdge(0));
    EXPECT_TRUE(hGraph.hasEdge(1));
    EXPECT_TRUE(hGraph.hasEdge(2));
    EXPECT_TRUE(hGraph.hasEdge(3));
    EXPECT_TRUE(hGraph.hasEdge(4));
    EXPECT_FALSE(hGraph.hasEdge(5));
    EXPECT_FALSE(hGraph.hasEdge(6));

    hGraph.removeEdge(0);
    hGraph.removeEdge(2);
    hGraph.addEdge();

    EXPECT_FALSE(hGraph.hasEdge(0));
    EXPECT_TRUE(hGraph.hasEdge(1));
    EXPECT_FALSE(hGraph.hasEdge(2));
    EXPECT_TRUE(hGraph.hasEdge(3));
    EXPECT_TRUE(hGraph.hasEdge(4));
    EXPECT_TRUE(hGraph.hasEdge(5));
    EXPECT_FALSE(hGraph.hasEdge(6));
}

TEST_F(HypergraphGTest, testRemoveEdge) {
    auto testGraph = [&](Hypergraph &hGraph) {
        count numEdges = hGraph.numberOfEdges();
        count maxEdgeId = numEdges;
        for (edgeid eid = 0; eid < maxEdgeId; ++eid) {
            hGraph.removeEdge(eid);
            --numEdges;
            EXPECT_EQ(hGraph.numberOfEdges(), numEdges);
            hGraph.forEdges(
                [&](edgeid edgeTest) { EXPECT_EQ(hGraph.hasEdge(edgeTest), edgeTest > eid); });
        }
    };

    Hypergraph hGraph = SimpleHypergraphGenerator(100, 100, 5).generate();
    testGraph(hGraph);
}

/** EDGE PROPERTIES **/

TEST_P(HypergraphGTest, testGetSetEdgeWeight) {
    Hypergraph hGraph = Hypergraph(0, 3, isWeightedHypergraph());

    EXPECT_DOUBLE_EQ(hGraph.getEdgeWeight(0), defaultEdgeWeight);

    hGraph.setEdgeWeight(1, 2.0);

    EXPECT_DOUBLE_EQ(hGraph.getEdgeWeight(1), isWeightedHypergraph() ? 2.0 : defaultEdgeWeight);

    hGraph.setEdgeWeight(2, -2.0);

    EXPECT_DOUBLE_EQ(hGraph.getEdgeWeight(2), isWeightedHypergraph() ? -2.0 : defaultEdgeWeight);
}

TEST_F(HypergraphGTest, testOrder) {
    Hypergraph hGraph = SimpleHypergraphGenerator(100, 100, 10, true).generate();
    hGraph.forEdges([&](edgeid eid) { EXPECT_EQ(10u, hGraph.order(eid)); });
}

/** EDGE ITERATORS **/

TEST_P(HypergraphGTest, testForEdges) {
    double epsilon = 1e-6;

    Hypergraph hGraph = Hypergraph(5, 4, true);
    hGraph.setEdgeWeight(0, 0.1);
    hGraph.setEdgeWeight(1, 0.2);
    hGraph.setEdgeWeight(2, 0.3);
    hGraph.setEdgeWeight(3, 0.4);

    hGraph.addNodeTo({0}, 0);
    hGraph.addNodeTo({0, 1}, 1);
    hGraph.addNodeTo({1, 2}, 2);
    hGraph.addNodeTo({2, 3}, 3);
    hGraph.addNodeTo({3}, 4);

    std::vector<bool> edgesSeen(4, false);
    edgeweight weightSum = 0;

    hGraph.forEdges([&](edgeid eid, edgeweight ew) {
        EXPECT_TRUE(hGraph.hasNode(static_cast<node>(eid), eid));
        EXPECT_TRUE(hGraph.hasNode(static_cast<node>(eid + 1), eid));
        EXPECT_EQ(hGraph.getEdgeWeight(eid), ew);

        edgesSeen[static_cast<index>(eid)] = true;

        if (hGraph.isWeighted()) {
            EXPECT_NEAR((eid + 1) * 0.1, ew, epsilon);
        } else {
            EXPECT_EQ(defaultEdgeWeight, ew);
        }
        weightSum += ew;
    });

    for (auto b : edgesSeen) {
        EXPECT_TRUE(b);
    }
    if (hGraph.isWeighted()) {
        EXPECT_NEAR(1.0, weightSum, epsilon);
    } else {
        EXPECT_NEAR(4 * defaultEdgeWeight, weightSum, epsilon);
    }
}

TEST_P(HypergraphGTest, testParallelForEdges) {
    double epsilon = 1e-6;

    Hypergraph hGraph = Hypergraph(5, 4, true);
    hGraph.setEdgeWeight(0, 0.1);
    hGraph.setEdgeWeight(1, 0.2);
    hGraph.setEdgeWeight(2, 0.3);
    hGraph.setEdgeWeight(3, 0.4);

    hGraph.addNodeTo({0}, 0);
    hGraph.addNodeTo({0, 1}, 1);
    hGraph.addNodeTo({1, 2}, 2);
    hGraph.addNodeTo({2, 3}, 3);
    hGraph.addNodeTo({3}, 4);

    std::vector<bool> edgesSeen(4, false);
    edgeweight weightSum = 0;

    hGraph.parallelForEdges([&](edgeid eid, edgeweight ew) {
        EXPECT_TRUE(hGraph.hasNode(static_cast<node>(eid), eid));
        EXPECT_TRUE(hGraph.hasNode(static_cast<node>(eid + 1), eid));
        EXPECT_EQ(hGraph.getEdgeWeight(eid), ew);

        edgesSeen[static_cast<index>(eid)] = true;

        if (hGraph.isWeighted()) {
            EXPECT_NEAR((eid + 1) * 0.1, ew, epsilon);
        } else {
            EXPECT_EQ(defaultEdgeWeight, ew);
        }
#pragma omp atomic
        weightSum += ew;
    });

    for (auto b : edgesSeen) {
        EXPECT_TRUE(b);
    }
    if (hGraph.isWeighted()) {
        EXPECT_NEAR(1.0, weightSum, epsilon);
    } else {
        EXPECT_NEAR(4 * defaultEdgeWeight, weightSum, epsilon);
    }
}

/* NEIGHBOR ITERATORS */

TEST_P(HypergraphGTest, testForNeighborsOf) {
    std::vector<node> visited;
    Hypergraph hGraph = Hypergraph(5, 5);
    hGraph.addNodesTo({0, 1}, 0);
    hGraph.addNodesTo({1}, 1);
    hGraph.addNodesTo({1, 2}, 2);
    hGraph.addNodesTo({1, 4}, 4);

    hGraph.forNeighborsOf(1, [&](node u) { visited.push_back(u); });

    EXPECT_EQ(3u, visited.size());
    EXPECT_EQ(0u, visited[0]);
    EXPECT_EQ(2u, visited[1]);
    EXPECT_EQ(4u, visited[2]);
}

} // namespace NetworKit
