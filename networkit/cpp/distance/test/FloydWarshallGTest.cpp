/*  FloydWarshallTestFixture.cpp
 *
 *	Created on: 19.02.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <networkit/distance/FloydWarshall.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace {
using ::testing::ElementsAre;
using ::testing::IsEmpty;

using GraphTypes = ::testing::Types<DynamicGraph<uint64_t, double>, DynamicGraph<int32_t, float>>;

template <class GraphType>
class FloydWarshallTestFixture : public testing::Test {
    using NodeType = typename GraphType::node_type;
    using EdgeWeightType = typename GraphType::edge_weight_type;

public:
    static constexpr EdgeWeightType maxDistance = std::numeric_limits<EdgeWeightType>::max();

    GraphType completeGraphK3() {
        GraphType graph(3, true);
        graph.addEdge(0, 1, 1);
        graph.addEdge(1, 2, 2);
        graph.addEdge(0, 2, 4);
        return graph;
    }

    GraphType undirectedGraphWithNegativeEdge() {
        GraphType graph(3, true);
        graph.addEdge(0, 1, 1);
        graph.addEdge(1, 2, 2);
        graph.addEdge(0, 2, -0.5);
        return graph;
    }

    GraphType directedGraphNegativeEdge() {
        GraphType graph(3, true, true);
        graph.addEdge(0, 1, 1);
        graph.addEdge(1, 2, -2);
        graph.addEdge(0, 2, 4);
        return graph;
    }

    GraphType completGraphK5NegativeEdges() {
        GraphType graph(5, true);
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

    GraphType disconnectedGraph() {
        GraphType graph(4, true);
        graph.addEdge(0, 1, 3);
        graph.addEdge(1, 2, 2);

        return graph;
    }

    GraphType directedGraphWithNegativeSelfLoop() {
        GraphType graph(5, true, true);
        graph.addEdge(0, 1, 3);
        graph.addEdge(1, 1, -2); // self-loop with negative cycle
        graph.addEdge(1, 2, 2);
        graph.addEdge(2, 3, 1);
        graph.addEdge(3, 4, 4);
        graph.addEdge(4, 0, 1);
        return graph;
    }

    void compareDistances(const std::vector<std::vector<EdgeWeightType>> &expectedDistances,
                          const FloydWarshall<GraphType> &testObject, const GraphType &graph) {
        for (NodeType source : graph.nodeRange()) {
            for (NodeType target : graph.nodeRange()) {
                EXPECT_EQ(testObject.getDistance(source, target), expectedDistances[source][target])
                    << "source = " << source << ", target = " << target;
                ;
            }
        }
    }

    void compareNodesOnShortestPaths(
        const std::vector<std::vector<std::vector<NodeType>>> &expectedPaths,
        const FloydWarshall<GraphType> &testObject, const GraphType &graph) {
        for (NodeType source : graph.nodeRange()) {
            for (NodeType target : graph.nodeRange()) {
                EXPECT_EQ(testObject.getNodesOnShortestPath(source, target),
                          expectedPaths[source][target])
                    << "source = " << source << ", target = " << target;
            }
        }
    }
};

TYPED_TEST_SUITE(FloydWarshallTestFixture, GraphTypes, /*comma required for variadic macro*/);

TYPED_TEST(FloydWarshallTestFixture, testConstructorThrowsUnweightedGraph) {
    const TypeParam graph(1, false);
    try {
        FloydWarshall<TypeParam> test(graph);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "The input graph is unweighted!");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TYPED_TEST(FloydWarshallTestFixture, testGetDistanceThrows) {
    const TypeParam graph(1, true);
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

TYPED_TEST(FloydWarshallTestFixture, testIsNodeInNegativeCycleThrows) {
    const TypeParam graph(1, true);
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

TYPED_TEST(FloydWarshallTestFixture, testGetNodesOnShortestPathThrows) {
    const TypeParam graph(2, true);
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

TYPED_TEST(FloydWarshallTestFixture, testGetDistanceCompleteGraphK3) {
    const auto graph = this->completeGraphK3();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<typename TypeParam::edge_weight_type>> expectedDistances{
        {0, 1, 3}, {1, 0, 2}, {3, 2, 0}};
    this->compareDistances(expectedDistances, test, graph);
}

TYPED_TEST(FloydWarshallTestFixture, testIsNodeInNegativeCycleCompleteGraphK3) {
    const auto graph = this->completeGraphK3();
    FloydWarshall test(graph);
    test.run();
    graph.forNodes([&test](const auto u) { EXPECT_FALSE(test.isNodeInNegativeCycle(u)); });
}

TYPED_TEST(FloydWarshallTestFixture, getNodesOnShortestPathCompleteGraphK3) {
    const auto graph = this->completeGraphK3();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<std::vector<typename TypeParam::node_type>>>
        expectedNodesOnShortestPaths{
            {{0}, {0, 1}, {0, 1, 2}}, {{1, 0}, {1}, {1, 2}}, {{2, 1, 0}, {2, 1}, {2}}};
    this->compareNodesOnShortestPaths(expectedNodesOnShortestPaths, test, graph);
}

TYPED_TEST(FloydWarshallTestFixture, testGetDistanceUndirectedGraphWithNegativeEdge) {
    // Undirected graph with one negative edge results in negative cycles for all edges
    const auto graph = this->undirectedGraphWithNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    using EdgeWeightType = typename TypeParam::edge_weight_type;
    constexpr EdgeWeightType expectedDistance{-std::numeric_limits<EdgeWeightType>::infinity()};
    for (const auto source : graph.nodeRange()) {
        for (const auto target : graph.nodeRange()) {
            EXPECT_EQ(test.getDistance(source, target), expectedDistance);
        }
    }
}

TYPED_TEST(FloydWarshallTestFixture, testIsNodeInNegativeCycleUndirectedGraphWithNegativeEdge) {
    // Undirected graph with one negative edge results in negative cycles for all edges
    const auto graph = this->undirectedGraphWithNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    for (const auto source : graph.nodeRange()) {
        EXPECT_TRUE(test.isNodeInNegativeCycle(source));
    }
}

TYPED_TEST(FloydWarshallTestFixture, getNodesOnShortestPathUndirectedGraphWithNegativeEdge) {
    // Undirected graph with one negative edge results in negative cycles for all edges
    const auto graph = this->undirectedGraphWithNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    for (const auto source : graph.nodeRange()) {
        for (const auto target : graph.nodeRange()) {
            EXPECT_THAT(test.getNodesOnShortestPath(source, target), IsEmpty());
        }
    }
}

TYPED_TEST(FloydWarshallTestFixture, testGetDistanceCompleteGraphK3NegativeEdge) {
    const auto graph = this->directedGraphNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<typename TypeParam::edge_weight_type>> expectedDistances{
        {0, 1, -1}, {this->maxDistance, 0, -2}, {this->maxDistance, this->maxDistance, 0}};
    this->compareDistances(expectedDistances, test, graph);
}

TYPED_TEST(FloydWarshallTestFixture, testIsNodeInNegativeEdgeCompleteGraphK3NegativeEdge) {
    const auto graph = this->directedGraphNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    for (const auto source : graph.nodeRange()) {
        EXPECT_FALSE(test.isNodeInNegativeCycle(source));
    }
}

TYPED_TEST(FloydWarshallTestFixture, testGetNodesOnShortestPathCompleteGraphK3NegativeEdge) {
    const auto graph = this->directedGraphNegativeEdge();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<std::vector<typename TypeParam::node_type>>>
        expectedNodesOnShortestPaths{{{0}, {0, 1}, {0, 1, 2}}, {{}, {1}, {1, 2}}, {{}, {}, {2}}};
    this->compareNodesOnShortestPaths(expectedNodesOnShortestPaths, test, graph);
}

TYPED_TEST(FloydWarshallTestFixture, testGetDistanceDisconnectedGraph) {
    const auto graph = this->disconnectedGraph();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<typename TypeParam::edge_weight_type>> expectedDistances{
        {0, 3, 5, this->maxDistance},
        {3, 0, 2, this->maxDistance},
        {5, 2, 0, this->maxDistance},
        {this->maxDistance, this->maxDistance, this->maxDistance, 0}};
    this->compareDistances(expectedDistances, test, graph);
}

TYPED_TEST(FloydWarshallTestFixture, testGetNodesOnShortestPathDisconnectedGraph) {
    const auto graph = this->disconnectedGraph();
    FloydWarshall test(graph);
    test.run();
    const std::vector<std::vector<std::vector<typename TypeParam::node_type>>>
        expectedNodesOnShortestPaths{{{0}, {0, 1}, {0, 1, 2}, {}},
                                     {{1, 0}, {1}, {1, 2}, {}},
                                     {{2, 1, 0}, {2, 1}, {2}, {}},
                                     {{}, {}, {}, {3}}};
    this->compareNodesOnShortestPaths(expectedNodesOnShortestPaths, test, graph);
}

TYPED_TEST(FloydWarshallTestFixture, testIsNodeInNegativeCycleDisconnectedGraph) {
    const auto graph = this->disconnectedGraph();
    FloydWarshall test(graph);
    test.run();
    for (const auto source : graph.nodeRange()) {
        EXPECT_FALSE(test.isNodeInNegativeCycle(source));
    }
}

TYPED_TEST(FloydWarshallTestFixture, testGetDistanceDirectedGraphWithNegativeSelfLoop) {
    const auto graph = this->directedGraphWithNegativeSelfLoop();
    FloydWarshall test(graph);
    test.run();
    using EdgeWeightType = typename TypeParam::edge_weight_type;
    constexpr EdgeWeightType expectedDistance{-std::numeric_limits<EdgeWeightType>::infinity()};
    for (const auto source : graph.nodeRange()) {
        for (const auto target : graph.nodeRange()) {
            EXPECT_EQ(test.getDistance(source, target), expectedDistance);
        }
    }
}

TYPED_TEST(FloydWarshallTestFixture, testIsNodeInNegativeCycleDirectedGraphWithNegativeSelfLoop) {
    const auto graph = this->directedGraphWithNegativeSelfLoop();
    FloydWarshall test(graph);
    test.run();
    for (const auto source : graph.nodeRange()) {
        EXPECT_TRUE(test.isNodeInNegativeCycle(source));
    }
}

TYPED_TEST(FloydWarshallTestFixture, testMultipleShortestDistancePaths) {
    TypeParam graph(11, true);
    // Shortest path, first case [0,10] with 5 nodes (inclusive)
    graph.addEdge(0, 1, 1);
    graph.addEdge(1, 2, 1);
    graph.addEdge(2, 3, 1);
    graph.addEdge(3, 10, 2);
    // Shortest path, second case [0,10] with 4 nodes (inclusive)
    graph.addEdge(0, 4, 1);
    graph.addEdge(4, 5, 1);
    graph.addEdge(5, 10, 3);
    // Shortest path, third case [0,10] with 6 nodes (inclusive)
    graph.addEdge(0, 6, 1);
    graph.addEdge(6, 7, 1);
    graph.addEdge(7, 8, 1);
    graph.addEdge(8, 9, 1);
    graph.addEdge(9, 10, 1);
    FloydWarshall test(graph);
    test.run();
    EXPECT_EQ(test.getDistance(0, 10), 5.0);
    EXPECT_THAT(test.getNodesOnShortestPath(0, 10), ElementsAre(0, 4, 5, 10));
}
} // namespace
} // namespace NetworKit
