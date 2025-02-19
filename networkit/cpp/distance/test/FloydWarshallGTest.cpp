/*  FloydWarshallGTest.cpp
 *
 *	Created on: 19.02.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <gtest/gtest.h>
#include <networkit/distance/FloydWarshall.hpp>

namespace NetworKit {

TEST(FloydWarshallTest, testConstructorThrowsUnweightedGraph) {
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

TEST(FloydWarshallTest, testConstructorThrowsInvalidDensityThreshold) {
    Graph graph(1, true);
    try {
        FloydWarshall test(graph, 0.0);
        FAIL() << "Expected std::invalid_argument";
    } catch (const std::invalid_argument &e) {
        EXPECT_STREQ(e.what(), "Invalid density threshold. Must be in range (0, 1].");
    } catch (...) {
        FAIL() << "Expected std::invalid_argument but got a different exception.";
    }
    try {
        FloydWarshall test(graph, 1.1);
        FAIL() << "Expected std::invalid_argument";
    } catch (const std::invalid_argument &e) {
        EXPECT_STREQ(e.what(), "Invalid density threshold. Must be in range (0, 1].");
    } catch (...) {
        FAIL() << "Expected std::invalid_argument but got a different exception.";
    }
}

TEST(FloydWarshallTest, testConstructorThrowsInvalidMaximumNumberOfNodes) {
    Graph graph(1, true);
    try {
        FloydWarshall test(graph, 0.5, 0);
        FAIL() << "Expected std::invalid_argument";
    } catch (const std::invalid_argument &e) {
        EXPECT_STREQ(e.what(), "Invalid maximum node count. Must be at least 1.");
    } catch (...) {
        FAIL() << "Expected std::invalid_argument but got a different exception.";
    }
}

TEST(FloydWarshallTest, testConstructorThrowsInvalidDensity) {
    Graph graph(2, true);
    constexpr double densityThreshold = 0.5;
    const std::string exceptionString = "Graph density is below user-defined density-threshold of: "
                                        + std::to_string(densityThreshold);
    try {
        FloydWarshall test(graph, densityThreshold);
        FAIL() << "Expected std::domain_error";
    } catch (const std::domain_error &e) {
        EXPECT_STREQ(e.what(), exceptionString.c_str());
    } catch (...) {
        FAIL() << "Expected std::domain_error but got a different exception.";
    }
}

TEST(FloydWarshallTest, testConstructorThrowsInvalidNumberOfNodes) {
    Graph graph(2, true);
    graph.addEdge(0, 0, 1);
    constexpr node maximumNumberOfNodes = 1;
    const std::string exceptionString = "Graph size exceeds user-defined max-node-limit of: "
                                        + std::to_string(maximumNumberOfNodes);
    try {
        constexpr double densityThreshold = 0.5;
        FloydWarshall test(graph, densityThreshold, maximumNumberOfNodes);
        FAIL() << "Expected std::domain_error";
    } catch (const std::domain_error &e) {
        EXPECT_STREQ(e.what(), exceptionString.c_str());
    } catch (...) {
        FAIL() << "Expected std::domain_error but got a different exception.";
    }
}

TEST(FloydWarshallTest, testGetDistanceThrows) {
    constexpr index numberOfNodes{3};
    Graph graph(numberOfNodes, true);
    graph.addEdge(0, 1, 1);
    graph.addEdge(1, 2, 2);
    graph.addEdge(0, 2, 4);
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

TEST(FloydWarshallTest, testGetDistance) {
    constexpr index numberOfNodes{3};
    Graph graph(numberOfNodes, true);
    graph.addEdge(0, 1, 1);
    graph.addEdge(1, 2, 2);
    graph.addEdge(0, 2, 4);
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<edgeweight>> expectedDistances{{0, 1, 3}, {1, 0, 2}, {3, 2, 0}};
    for (node u = 0; u < numberOfNodes; ++u) {
        for (node v = 0; v < numberOfNodes; ++v) {
            EXPECT_EQ(test.getDistance(u, v), expectedDistances[u][v]) << "u = " << u << ", v = " << v;
        }
    }
}

TEST(FloydWarshallTest, testGetAllDistancesThrows) {
    constexpr index numberOfNodes{3};
    Graph graph(numberOfNodes, true);
    graph.addEdge(0, 1, 1);
    graph.addEdge(1, 2, 2);
    graph.addEdge(0, 2, 4);
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

TEST(FloydWarshallTest, testAllGetDistances) {
    constexpr index numberOfNodes{3};
    Graph graph(numberOfNodes, true);
    graph.addEdge(0, 1, 1);
    graph.addEdge(1, 2, 2);
    graph.addEdge(0, 2, 4);
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<edgeweight>> expectedDistances{{0, 1, 3}, {1, 0, 2}, {3, 2, 0}};
    const auto distances = test.getAllDistances();
    for (node u = 0; u < numberOfNodes; ++u) {
        for (node v = 0; v < numberOfNodes; ++v) {
            EXPECT_EQ(distances[u][v], expectedDistances[u][v]) << "u = " << u << ", v = " << v;
        }
    }
}


} // namespace NetworKit
