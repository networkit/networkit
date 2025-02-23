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

    Graph completeGraphK3NegativeCycle() {
        Graph graph(3, true);
        graph.addEdge(0, 1, 1);
        graph.addEdge(1, 2, 2);
        graph.addEdge(0, 2, -4);
        return graph;
    }

    Graph directedGraphNegativeEdge() {
        Graph graph(3, true, true);
        graph.addEdge(0, 1, 1);
        graph.addEdge(1, 2, -2);
        graph.addEdge(0, 2, 4);
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

TEST_F(FloydWarshallGTest, testGetAllDistancesThrows) {
    Graph graph(1, true);
    FloydWarshall test(graph);
    try {
        test.getAllDistances();
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

TEST_F(FloydWarshallGTest, testAllGetDistancesCompleteGraphK3) {
    auto graph = completeGraphK3();
    FloydWarshall test(graph);    constexpr edgeweight maxDistance = std::numeric_limits<edgeweight>::max();

    test.run();
    const std::vector<std::vector<edgeweight>> expectedDistances{{0, 1, 3}, {1, 0, 2}, {3, 2, 0}};
    const auto & distances = test.getAllDistances();
    EXPECT_EQ(distances, expectedDistances);
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
            EXPECT_EQ(test.getNodesOnShortestPath(source, target), expectedNodesOnShortestPaths[source][target]);
        }
    }
}

TEST_F(FloydWarshallGTest, testGetDistanceCompleteGraphK3NegativeCycle) {
    auto graph = completeGraphK3NegativeCycle();
    FloydWarshall test(graph);
    test.run();
    constexpr edgeweight expectedDistance{-std::numeric_limits<edgeweight>::infinity()};
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        for (node target = 0; target < graph.numberOfNodes(); ++target) {
            EXPECT_EQ(test.getDistance(source, target), expectedDistance);
        }
    }
}

TEST_F(FloydWarshallGTest, testGetAllDistancesCompleteGraphK3NegativeCycle) {
    auto graph = completeGraphK3NegativeCycle();
    FloydWarshall test(graph);
    test.run();
    constexpr edgeweight expectedDistance{-std::numeric_limits<edgeweight>::infinity()};
    const auto & distances = test.getAllDistances();
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        for (node target = 0; target < graph.numberOfNodes(); ++target) {
            EXPECT_EQ(distances[source][target], expectedDistance);
        }
    }
}

TEST_F(FloydWarshallGTest, testIsNodeInNegativeCycleCompleteGraphK3NegativeCycle) {
    auto graph = completeGraphK3NegativeCycle();
    FloydWarshall test(graph);
    test.run();
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        EXPECT_TRUE(test.isNodeInNegativeCycle(source));
    }
}

TEST_F(FloydWarshallGTest, getNodesOnShortestPathCompleteGraphK3NegativeCycle) {
    auto graph = completeGraphK3NegativeCycle();
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
    const std::vector<std::vector<edgeweight>> expectedDistances{{0, 1, -1}, {maxDistance, 0, -2}, {maxDistance,maxDistance, 0}};
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        for (node target = 0; target < graph.numberOfNodes(); ++target) {
            EXPECT_EQ(test.getDistance(source, target), expectedDistances[source][target]) << "source " << source << " target " << target;
        }
    }
}

TEST_F(FloydWarshallGTest, testGetAllDistancesCompleteGraphK3NegativeEdge) {
    auto graph = directedGraphNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<edgeweight>> expectedDistances{{0, 1, -1}, {maxDistance, 0, -2}, {maxDistance,maxDistance, 0}};
    const auto & distances = test.getAllDistances();
    EXPECT_EQ(distances, expectedDistances);
}

TEST_F(FloydWarshallGTest, testIsNodeInNegativeEdgeCompleteGraphK3NegativeEdge) {
    auto graph = directedGraphNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        EXPECT_FALSE(test.isNodeInNegativeCycle(source));
    }
}

TEST_F(FloydWarshallGTest, getNodesOnShortestPathCompleteGraphK3NegativeEdge) {
    auto graph = directedGraphNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<std::vector<node>>> expectedNodesOnShortestPaths{
            {{0}, {0, 1}, {0, 1, 2}}, {{}, {1}, {1, 2}}, {{}, {}, {2}}};
    for (node source = 0; source < graph.numberOfNodes(); ++source) {
        for (node target = 0; target < graph.numberOfNodes(); ++target) {
            EXPECT_EQ(test.getNodesOnShortestPath(source, target), expectedNodesOnShortestPaths[source][target]) << source << ", " << target;
        }
    }
}

TEST_F(FloydWarshallGTest, testDisconnectedGraph) {
    Graph graph(4, true);
    graph.addEdge(0, 1, 3);
    graph.addEdge(1, 2, 2);
    // No edges to/from node 3
    FloydWarshall test(graph);
    test.run();
    EXPECT_EQ(test.getDistance(0, 3), std::numeric_limits<edgeweight>::max());
    EXPECT_EQ(test.getDistance(1, 3), std::numeric_limits<edgeweight>::max());
    EXPECT_EQ(test.getDistance(2, 3), std::numeric_limits<edgeweight>::max());
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

TEST_F(FloydWarshallGTest, testSelfLoops) {
    Graph graph(3, true);
    graph.addEdge(0, 0, -2);
    graph.addEdge(1, 1, 3);
    FloydWarshall test(graph);
    test.run();
    EXPECT_EQ(test.getDistance(0, 0), -std::numeric_limits<edgeweight>::infinity());
    EXPECT_EQ(test.getDistance(1, 1), 3);
}

} // namespace NetworKit
