/*  FloydWarshallGTest.cpp
 *
 *	Created on: 19.02.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <gtest/gtest.h>
#include <networkit/distance/FloydWarshall.hpp>

namespace NetworKit {

class FloydWarshallGTest : public testing::Test {
public:
    static constexpr edgeweight maxDistance = std::numeric_limits<edgeweight>::max();
    Graph completeGraphK3() {
        Graph graph(3, true);
        graph.addEdge(0, 1, 1);
        graph.addEdge(1, 2, 2);
        graph.addEdge(0, 2, 4);
        return graph;
    }

    Graph undirectedGraphWithNegativeEdge() {
        Graph graph(3, true);
        graph.addEdge(0, 1, 1);
        graph.addEdge(1, 2, 2);
        graph.addEdge(0, 2, -0.5);
        return graph;
    }

    Graph directedGraphNegativeEdge() {
        Graph graph(3, true, true);
        graph.addEdge(0, 1, 1);
        graph.addEdge(1, 2, -2);
        graph.addEdge(0, 2, 4);
        return graph;
    }

    Graph completGraphK5NegativeEdges() {
        Graph graph(5, true);
        graph.addEdge(0, 1, -3);
        graph.addEdge(2, 3, -2);
        graph.addEdge(0, 2, 5);
        graph.addEdge(0, 3, 4);
        graph.addEdge(0, 4, 6);
        graph.addEdge(1, 2, 2);
        graph.addEdge(1, 3, 3);
        graph.addEdge(1, 4, 7);
        graph.addEdge(2, 4, 4);
        graph.addEdge(3, 4, 5);

        return graph;
    }

    Graph disconnectedGraph() {
        Graph graph(4, true);
        graph.addEdge(0, 1, 3);
        graph.addEdge(1, 2, 2);

        return graph;
    }
    //        ___
    //       |   | self-loop on node 1
    // 0---->1<---
    // A     |
    // |     |  --->5
    // |     V /    |
    // |     2      |
    // |     | \    |
    // |     V   \  |
    // 4<----3    \ V
    //              6
    Graph directedGraphWithNegativeSelfLoop() {
        Graph graph(7, true, true);
        graph.addEdge(0, 1, 3);
        graph.addEdge(1, 1, -2); // self-loop with negative cycle
        graph.addEdge(1, 2, 2);
        graph.addEdge(2, 3, 1);
        graph.addEdge(3, 4, 4);
        graph.addEdge(4, 0, 1);
        // Graph without self-loop
        graph.addEdge(2, 5, 2);
        graph.addEdge(5, 6, 3);
        graph.addEdge(6, 2, 3);

        return graph;
    }

};

TEST_F(FloydWarshallGTest, testConstructorThrowsUnweightedGraph) {
    Graph graph(1, false);
    try {
        FloydWarshall test(graph);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "The input graph is unweighted!");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(FloydWarshallGTest, testGetDistanceThrows) {
    Graph graph(1, true);
    FloydWarshall test(graph);
    try {
        test.getDistance(0, 1);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Error, run must be called first");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(FloydWarshallGTest, testIsNodeInNegativeCycleThrows) {
    Graph graph(1, true);
    FloydWarshall test(graph);
    try {
        test.isNodeInNegativeCycle(0);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Error, run must be called first");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(FloydWarshallGTest, testGetNodesOnShortestPathThrows) {
    Graph graph(2, true);
    FloydWarshall test(graph);
    try {
        test.getNodesOnShortestPath(0, 1);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Error, run must be called first");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(FloydWarshallGTest, testGetDistanceCompleteGraphK3) {
    auto graph = completeGraphK3();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<edgeweight>> expectedDistances{{0, 1, 3}, {1, 0, 2}, {3, 2, 0}};
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        for (node target = 0; target < graph.numberOfNodes(); ++target) {
            EXPECT_EQ(test.getDistance(source, target), expectedDistances[source][target])
                << "source = " << source << ", target = " << target;
        }
    }
}

TEST_F(FloydWarshallGTest, testIsNodeInNegativeCycleCompleteGraphK3) {
    auto graph = completeGraphK3();
    FloydWarshall test(graph);
    test.run();
    for (node u = 0; u < graph.numberOfNodes(); ++u) {
        EXPECT_FALSE(test.isNodeInNegativeCycle(u));
    }
}

TEST_F(FloydWarshallGTest, getNodesOnShortestPathCompleteGraphK3) {
    auto graph = completeGraphK3();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<std::vector<node>>> expectedNodesOnShortestPaths{
        {{0}, {0, 1}, {0, 1, 2}}, {{1, 0}, {1}, {1, 2}}, {{2, 1, 0}, {2, 1}, {2}}};
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        for (node target = 0; target < graph.numberOfNodes(); ++target) {
            EXPECT_EQ(test.getNodesOnShortestPath(source, target),
                      expectedNodesOnShortestPaths[source][target]);
        }
    }
}

TEST_F(FloydWarshallGTest, testGetDistanceUndirectedGraphWithNegativeEdge) {
    // Undirected graph with one negative edge results in negative cycles for all edges
    auto graph = undirectedGraphWithNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    constexpr edgeweight expectedDistance{-std::numeric_limits<edgeweight>::infinity()};
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        for (node target = 0; target < graph.numberOfNodes(); ++target) {
            EXPECT_EQ(test.getDistance(source, target), expectedDistance);
        }
    }
}

TEST_F(FloydWarshallGTest, testIsNodeInNegativeCycleUndirectedGraphWithNegativeEdge) {
    // Undirected graph with one negative edge results in negative cycles for all edges
    auto graph = undirectedGraphWithNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        EXPECT_TRUE(test.isNodeInNegativeCycle(source));
    }
}

TEST_F(FloydWarshallGTest, getNodesOnShortestPathUndirectedGraphWithNegativeEdge) {
    // Undirected graph with one negative edge results in negative cycles for all edges
    auto graph = undirectedGraphWithNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        for (node target = 0; target < graph.numberOfNodes(); ++target) {
            EXPECT_EQ(test.getNodesOnShortestPath(source, target), std::vector<node>{});
        }
    }
}

TEST_F(FloydWarshallGTest, testGetDistanceCompleteGraphK3NegativeEdge) {
    auto graph = directedGraphNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<edgeweight>> expectedDistances{
        {0, 1, -1}, {maxDistance, 0, -2}, {maxDistance, maxDistance, 0}};
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        for (node target = 0; target < graph.numberOfNodes(); ++target) {
            EXPECT_EQ(test.getDistance(source, target), expectedDistances[source][target])
                << "source " << source << " target " << target;
        }
    }
}

TEST_F(FloydWarshallGTest, testIsNodeInNegativeEdgeCompleteGraphK3NegativeEdge) {
    auto graph = directedGraphNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        EXPECT_FALSE(test.isNodeInNegativeCycle(source));
    }
}

TEST_F(FloydWarshallGTest, testGetNodesOnShortestPathCompleteGraphK3NegativeEdge) {
    auto graph = directedGraphNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<std::vector<node>>> expectedNodesOnShortestPaths{
        {{0}, {0, 1}, {0, 1, 2}}, {{}, {1}, {1, 2}}, {{}, {}, {2}}};
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        for (node target = 0; target < graph.numberOfNodes(); ++target) {
            EXPECT_EQ(test.getNodesOnShortestPath(source, target),
                      expectedNodesOnShortestPaths[source][target])
                << source << ", " << target;
        }
    }
}

TEST_F(FloydWarshallGTest, testGetDistanceDisconnectedGraph) {
    auto graph = disconnectedGraph();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<edgeweight>> expectedDistances{
        {0, 3, 5, maxDistance},
        {3, 0, 2, maxDistance},
        {5, 2, 0, maxDistance},
        {maxDistance, maxDistance, maxDistance, 0}};
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        for (node target = 0; target < graph.numberOfNodes(); ++target) {
            EXPECT_EQ(test.getDistance(source, target), expectedDistances[source][target])
                << "source " << source << " target " << target;
        }
    }
}



TEST_F(FloydWarshallGTest, testGetNodesOnShortestPathDisconnectedGraph) {
    auto graph = disconnectedGraph();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<std::vector<node>>> expectedNodesOnShortestPaths{
            {{0}, {0, 1}, {0, 1, 2}, {}}, {{1, 0}, {1}, {1, 2}, {}}, {{2, 1, 0}, {2, 1}, {2}, {}}
                , {{}, {}, {}, {3}}};
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        for (node target = 0; target < graph.numberOfNodes(); ++target) {
            EXPECT_EQ(test.getNodesOnShortestPath(source, target),
                      expectedNodesOnShortestPaths[source][target])
                << source << ", " << target;
        }
    }
}

TEST_F(FloydWarshallGTest, testIsNodeInNegativeCycleDisconnectedGraph) {
    auto graph = disconnectedGraph();
    FloydWarshall test(graph);
    test.run();
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        EXPECT_FALSE(test.isNodeInNegativeCycle(source));
    }
}



TEST_F(FloydWarshallGTest, testZeroWeightEdges) {
    Graph graph(3, true);
    graph.addEdge(0, 1, 0.0);
    graph.addEdge(1, 2, 0.0);
    graph.addEdge(0, 2, 3.0);
    FloydWarshall test(graph);
    test.run();
    EXPECT_EQ(test.getDistance(0, 1), 0.0);
    EXPECT_EQ(test.getDistance(1, 2), 0.0);
    EXPECT_EQ(test.getDistance(0, 2), 0.0);
}


} // namespace NetworKit
