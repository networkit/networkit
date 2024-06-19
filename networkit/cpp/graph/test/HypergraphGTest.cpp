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
    hGraph.forNodes([&](node u) { ASSERT_EQ(10u, hGraph.degree(u)); });
}

TEST_P(HypergraphGTest, testNeighborsWithoutDuplicates) {
    Hypergraph hGraph = Hypergraph(3, 3);
    hGraph.addNodesTo({0, 1}, 0);
    hGraph.addNodesTo({1}, 1);
    hGraph.addNodesTo({1, 2}, 2);

    auto neighbors = hGraph.getNeighbors(1);

    ASSERT_THAT(neighbors, testing::ElementsAre(0, 2));
}

TEST_P(HypergraphGTest, testNeighborsWithDuplicates) {
    Hypergraph hGraph = Hypergraph(3, 3);
    hGraph.addNodesTo({0, 1}, 0);
    hGraph.addNodesTo({1}, 1);
    hGraph.addNodesTo({0, 1}, 2);

    auto neighbors = hGraph.getNeighbors(1);

    ASSERT_EQ(neighbors.size(), 1);
    ASSERT_THAT(neighbors, testing::ElementsAre(0));
}

/** EDGE MODIFIERS **/

TEST_P(HypergraphGTest, testaddEdge) {
    Hypergraph hGraph = Hypergraph();

    ASSERT_FALSE(hGraph.hasEdge(0));
    ASSERT_FALSE(hGraph.hasEdge(1));
    ASSERT_EQ(0u, hGraph.numberOfEdges());

    hGraph.addEdge();
    ASSERT_TRUE(hGraph.hasEdge(0));
    ASSERT_FALSE(hGraph.hasEdge(1));
    ASSERT_EQ(1u, hGraph.numberOfEdges());

    Hypergraph hGraph2 = Hypergraph(0, 2);
    ASSERT_TRUE(hGraph2.hasEdge(0));
    ASSERT_TRUE(hGraph2.hasEdge(1));
    ASSERT_FALSE(hGraph2.hasEdge(2));
    ASSERT_EQ(2u, hGraph2.numberOfEdges());

    hGraph2.addEdge();
    hGraph2.addEdge();
    ASSERT_TRUE(hGraph2.hasEdge(2));
    ASSERT_TRUE(hGraph2.hasEdge(3));
    ASSERT_FALSE(hGraph2.hasEdge(4));
    ASSERT_EQ(4u, hGraph2.numberOfEdges());
}

TEST_P(HypergraphGTest, testaddEdgeWithNodes) {
    Hypergraph hGraph = Hypergraph(2);

    hGraph.addEdge({0, 1});
    ASSERT_TRUE(hGraph.hasEdge(0));
    ASSERT_EQ(1u, hGraph.numberOfEdges());
}

TEST_P(HypergraphGTest, testaddEdgeWithNodesAddMissing) {
    Hypergraph hGraph = Hypergraph();

    hGraph.addEdge({0, 1}, true);
    ASSERT_TRUE(hGraph.hasEdge(0));
    ASSERT_TRUE(hGraph.hasNode(0, 0));
    ASSERT_TRUE(hGraph.hasNode(1, 0));
    ASSERT_EQ(1u, hGraph.numberOfEdges());
}

/** EDGE ITERATORS **/

TEST_P(HypergraphGTest, testForEdges) {
    Hypergraph hGraph = Hypergraph(5, 4);
    hGraph.addNodeTo({0}, 0);
    hGraph.addNodeTo({0, 1}, 1);
    hGraph.addNodeTo({1, 2}, 2);
    hGraph.addNodeTo({2, 3}, 3);
    hGraph.addNodeTo({3}, 4);

    std::vector<bool> edgesSeen(4, false);

    hGraph.forEdges([&](edgeid eid) {
        ASSERT_TRUE(hGraph.hasNode(static_cast<node>(eid), eid));
        ASSERT_TRUE(hGraph.hasNode(static_cast<node>(eid + 1), eid));
        edgesSeen[static_cast<index>(eid)] = true;
    });

    for (auto b : edgesSeen) {
        ASSERT_TRUE(b);
    }
}

TEST_P(HypergraphGTest, testForWeightedEdges) {
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
        ASSERT_TRUE(hGraph.hasNode(static_cast<node>(eid), eid));
        ASSERT_TRUE(hGraph.hasNode(static_cast<node>(eid + 1), eid));
        ASSERT_EQ(hGraph.getEdgeWeight(eid), ew);

        edgesSeen[static_cast<index>(eid)] = true;

        if (hGraph.isWeighted()) {
            ASSERT_NEAR((eid + 1) * 0.1, ew, epsilon);
        } else {
            ASSERT_EQ(defaultEdgeWeight, ew);
        }
        weightSum += ew;
    });

    for (auto b : edgesSeen) {
        ASSERT_TRUE(b);
    }
    if (hGraph.isWeighted()) {
        ASSERT_NEAR(1.0, weightSum, epsilon);
    } else {
        ASSERT_NEAR(4 * defaultEdgeWeight, weightSum, epsilon);
    }
}

// TEST_P(GraphGTest, testParallelForWeightedEdges) {
//     count n = 4;
//     Graph G = createGraph(n);
//     G.forNodePairs([&](node u, node v) { G.addEdge(u, v, 1.0); });

//     edgeweight weightSum = 0.0;
//     G.parallelForEdges([&](node, node, edgeweight ew) {
// #pragma omp atomic
//         weightSum += ew;
//     });

//     ASSERT_EQ(6.0, weightSum) << "sum of edge weights should be 6 in every case";
// }

// TEST_P(GraphGTest, testParallelForEdges) {
//     count n = 4;
//     Graph G = createGraph(n);
//     G.forNodePairs([&](node u, node v) { G.addEdge(u, v); });

//     edgeweight weightSum = 0.0;
//     G.parallelForEdges([&](node, node) {
// #pragma omp atomic
//         weightSum += 1;
//     });

//     ASSERT_EQ(6.0, weightSum) << "sum of edge weights should be 6 in every case";
// }

/* NEIGHBOR ITERATORS */

TEST_P(HypergraphGTest, testForNeighborsOf) {
    std::vector<node> visited;
    Hypergraph hGraph = Hypergraph(5, 5);
    hGraph.addNodesTo({0, 1}, 0);
    hGraph.addNodesTo({1}, 1);
    hGraph.addNodesTo({1, 2}, 2);
    hGraph.addNodesTo({1, 4}, 4);

    hGraph.forNeighborsOf(1, [&](node u) { visited.push_back(u); });

    ASSERT_EQ(3u, visited.size());
    ASSERT_EQ(0u, visited[0]);
    ASSERT_EQ(2u, visited[1]);
    ASSERT_EQ(4u, visited[2]);
}

} // namespace NetworKit
