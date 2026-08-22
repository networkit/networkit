/*  FloydWarshallGTest.cpp
 *
 *	Created on: 19.02.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <limits>
#include <stdexcept>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/distance/FloydWarshall.hpp>
#include <networkit/graph/AdjListGraph.hpp>

namespace NetworKit {
namespace {

template <class NodeT_, class EdgeWeightT_>
struct FloydWarshallConfig {
    using NodeT = NodeT_;
    using EdgeWeightT = EdgeWeightT_;
};

template <class TestT>
class FloydWarshallGTest : public testing::Test {
public:
    using NodeT = typename TestT::NodeT;
    using EdgeWeightT = typename TestT::EdgeWeightT;
    using GraphT = AdjListGraph<NodeT, EdgeWeightT>;
    using FloydWarshallT = GenericFloydWarshall<GraphT>;

    static constexpr edgeweight maxDistance = std::numeric_limits<edgeweight>::max();

    GraphT completeGraphK3() const {
        GraphT graph(3, true);
        graph.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{1});
        graph.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{2});
        graph.addEdge(NodeT{0}, NodeT{2}, EdgeWeightT{4});
        return graph;
    }

    GraphT undirectedGraphWithNegativeEdge() const {
        GraphT graph(3, true);
        graph.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{1});
        graph.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{2});
        graph.addEdge(NodeT{0}, NodeT{2}, EdgeWeightT{-1});
        return graph;
    }

    GraphT directedGraphNegativeEdge() const {
        GraphT graph(3, true, true);
        graph.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{1});
        graph.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{-2});
        graph.addEdge(NodeT{0}, NodeT{2}, EdgeWeightT{4});
        return graph;
    }

    GraphT disconnectedGraph() const {
        GraphT graph(4, true);
        graph.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{3});
        graph.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{2});
        return graph;
    }

    GraphT directedGraphWithNegativeSelfLoop() const {
        GraphT graph(5, true, true);
        graph.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{3});
        graph.addEdge(NodeT{1}, NodeT{1}, EdgeWeightT{-2});
        graph.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{2});
        graph.addEdge(NodeT{2}, NodeT{3}, EdgeWeightT{1});
        graph.addEdge(NodeT{3}, NodeT{4}, EdgeWeightT{4});
        graph.addEdge(NodeT{4}, NodeT{0}, EdgeWeightT{1});
        return graph;
    }

    std::vector<std::vector<edgeweight>> distancesFromGetDistance(const FloydWarshallT &testObject,
                                                                  count numberOfNodes) const {
        std::vector<std::vector<edgeweight>> distances(numberOfNodes,
                                                       std::vector<edgeweight>(numberOfNodes));
        for (count source = 0; source < numberOfNodes; ++source) {
            for (count target = 0; target < numberOfNodes; ++target) {
                distances[source][target] =
                    testObject.getDistance(static_cast<NodeT>(source), static_cast<NodeT>(target));
            }
        }
        return distances;
    }

    std::vector<std::vector<std::vector<NodeT>>>
    pathsFromGetNodesOnShortestPath(const FloydWarshallT &testObject, count numberOfNodes) const {
        std::vector<std::vector<std::vector<NodeT>>> paths(
            numberOfNodes, std::vector<std::vector<NodeT>>(numberOfNodes));
        for (count source = 0; source < numberOfNodes; ++source) {
            for (count target = 0; target < numberOfNodes; ++target) {
                paths[source][target] = testObject.getNodesOnShortestPath(
                    static_cast<NodeT>(source), static_cast<NodeT>(target));
            }
        }
        return paths;
    }

    std::vector<bool> negativeCycleFlags(const FloydWarshallT &testObject,
                                         count numberOfNodes) const {
        std::vector<bool> flags(numberOfNodes);
        for (count u = 0; u < numberOfNodes; ++u) {
            flags[u] = testObject.isNodeInNegativeCycle(static_cast<NodeT>(u));
        }
        return flags;
    }

    void compareDistances(const std::vector<std::vector<edgeweight>> &expectedDistances,
                          const FloydWarshallT &testObject) const {
        EXPECT_THAT(distancesFromGetDistance(testObject, expectedDistances.size()),
                    testing::ElementsAreArray(expectedDistances));
    }

    void
    compareNodesOnShortestPaths(const std::vector<std::vector<std::vector<NodeT>>> &expectedPaths,
                                const FloydWarshallT &testObject) const {
        EXPECT_THAT(pathsFromGetNodesOnShortestPath(testObject, expectedPaths.size()),
                    testing::ElementsAreArray(expectedPaths));
    }
};

TYPED_TEST_SUITE_P(FloydWarshallGTest);

TYPED_TEST_P(FloydWarshallGTest, testConstructorThrowsUnweightedGraph) {
    using GraphT = typename TestFixture::GraphT;
    using FloydWarshallT = typename TestFixture::FloydWarshallT;

    GraphT graph(1, false);
    try {
        FloydWarshallT test(graph);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "The input graph is unweighted!");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TYPED_TEST_P(FloydWarshallGTest, testGetDistanceThrows) {
    using GraphT = typename TestFixture::GraphT;
    using FloydWarshallT = typename TestFixture::FloydWarshallT;
    using NodeT = typename TestFixture::NodeT;

    GraphT graph(1, true);
    FloydWarshallT test(graph);
    try {
        test.getDistance(NodeT{0}, NodeT{1});
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Error, run must be called first");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TYPED_TEST_P(FloydWarshallGTest, testIsNodeInNegativeCycleThrows) {
    using GraphT = typename TestFixture::GraphT;
    using FloydWarshallT = typename TestFixture::FloydWarshallT;
    using NodeT = typename TestFixture::NodeT;

    GraphT graph(1, true);
    FloydWarshallT test(graph);
    try {
        test.isNodeInNegativeCycle(NodeT{0});
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Error, run must be called first");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TYPED_TEST_P(FloydWarshallGTest, testGetNodesOnShortestPathThrows) {
    using GraphT = typename TestFixture::GraphT;
    using FloydWarshallT = typename TestFixture::FloydWarshallT;
    using NodeT = typename TestFixture::NodeT;

    GraphT graph(2, true);
    FloydWarshallT test(graph);
    try {
        test.getNodesOnShortestPath(NodeT{0}, NodeT{1});
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Error, run must be called first");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TYPED_TEST_P(FloydWarshallGTest, testGetDistanceCompleteGraphK3) {
    auto graph = this->completeGraphK3();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    const std::vector<std::vector<edgeweight>> expectedDistances{{0, 1, 3}, {1, 0, 2}, {3, 2, 0}};
    this->compareDistances(expectedDistances, test);
}

TYPED_TEST_P(FloydWarshallGTest, testIsNodeInNegativeCycleCompleteGraphK3) {
    auto graph = this->completeGraphK3();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    EXPECT_THAT(this->negativeCycleFlags(test, graph.numberOfNodes()), testing::Each(false));
}

TYPED_TEST_P(FloydWarshallGTest, getNodesOnShortestPathCompleteGraphK3) {
    using NodeT = typename TestFixture::NodeT;

    auto graph = this->completeGraphK3();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    const std::vector<std::vector<std::vector<NodeT>>> expectedNodesOnShortestPaths{
        {{NodeT{0}}, {NodeT{0}, NodeT{1}}, {NodeT{0}, NodeT{1}, NodeT{2}}},
        {{NodeT{1}, NodeT{0}}, {NodeT{1}}, {NodeT{1}, NodeT{2}}},
        {{NodeT{2}, NodeT{1}, NodeT{0}}, {NodeT{2}, NodeT{1}}, {NodeT{2}}}};
    this->compareNodesOnShortestPaths(expectedNodesOnShortestPaths, test);
}

TYPED_TEST_P(FloydWarshallGTest, testGetDistanceUndirectedGraphWithNegativeEdge) {
    auto graph = this->undirectedGraphWithNegativeEdge();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    constexpr edgeweight expectedDistance{-std::numeric_limits<edgeweight>::infinity()};
    EXPECT_THAT(this->distancesFromGetDistance(test, graph.numberOfNodes()),
                testing::Each(testing::Each(expectedDistance)));
}

TYPED_TEST_P(FloydWarshallGTest, testIsNodeInNegativeCycleUndirectedGraphWithNegativeEdge) {
    auto graph = this->undirectedGraphWithNegativeEdge();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    EXPECT_THAT(this->negativeCycleFlags(test, graph.numberOfNodes()), testing::Each(true));
}

TYPED_TEST_P(FloydWarshallGTest, getNodesOnShortestPathUndirectedGraphWithNegativeEdge) {
    auto graph = this->undirectedGraphWithNegativeEdge();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    EXPECT_THAT(this->pathsFromGetNodesOnShortestPath(test, graph.numberOfNodes()),
                testing::Each(testing::Each(testing::IsEmpty())));
}

TYPED_TEST_P(FloydWarshallGTest, testGetDistanceCompleteGraphK3NegativeEdge) {
    auto graph = this->directedGraphNegativeEdge();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    const std::vector<std::vector<edgeweight>> expectedDistances{
        {0, 1, -1}, {this->maxDistance, 0, -2}, {this->maxDistance, this->maxDistance, 0}};
    this->compareDistances(expectedDistances, test);
}

TYPED_TEST_P(FloydWarshallGTest, testIsNodeInNegativeEdgeCompleteGraphK3NegativeEdge) {
    auto graph = this->directedGraphNegativeEdge();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    EXPECT_THAT(this->negativeCycleFlags(test, graph.numberOfNodes()), testing::Each(false));
}

TYPED_TEST_P(FloydWarshallGTest, testGetNodesOnShortestPathCompleteGraphK3NegativeEdge) {
    using NodeT = typename TestFixture::NodeT;

    auto graph = this->directedGraphNegativeEdge();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    const std::vector<std::vector<std::vector<NodeT>>> expectedNodesOnShortestPaths{
        {{NodeT{0}}, {NodeT{0}, NodeT{1}}, {NodeT{0}, NodeT{1}, NodeT{2}}},
        {{}, {NodeT{1}}, {NodeT{1}, NodeT{2}}},
        {{}, {}, {NodeT{2}}}};
    this->compareNodesOnShortestPaths(expectedNodesOnShortestPaths, test);
}

TYPED_TEST_P(FloydWarshallGTest, testGetDistanceDisconnectedGraph) {
    auto graph = this->disconnectedGraph();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    const std::vector<std::vector<edgeweight>> expectedDistances{
        {0, 3, 5, this->maxDistance},
        {3, 0, 2, this->maxDistance},
        {5, 2, 0, this->maxDistance},
        {this->maxDistance, this->maxDistance, this->maxDistance, 0}};
    this->compareDistances(expectedDistances, test);
}

TYPED_TEST_P(FloydWarshallGTest, testGetNodesOnShortestPathDisconnectedGraph) {
    using NodeT = typename TestFixture::NodeT;

    auto graph = this->disconnectedGraph();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    const std::vector<std::vector<std::vector<NodeT>>> expectedNodesOnShortestPaths{
        {{NodeT{0}}, {NodeT{0}, NodeT{1}}, {NodeT{0}, NodeT{1}, NodeT{2}}, {}},
        {{NodeT{1}, NodeT{0}}, {NodeT{1}}, {NodeT{1}, NodeT{2}}, {}},
        {{NodeT{2}, NodeT{1}, NodeT{0}}, {NodeT{2}, NodeT{1}}, {NodeT{2}}, {}},
        {{}, {}, {}, {NodeT{3}}}};
    this->compareNodesOnShortestPaths(expectedNodesOnShortestPaths, test);
}

TYPED_TEST_P(FloydWarshallGTest, testIsNodeInNegativeCycleDisconnectedGraph) {
    auto graph = this->disconnectedGraph();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    EXPECT_THAT(this->negativeCycleFlags(test, graph.numberOfNodes()), testing::Each(false));
}

TYPED_TEST_P(FloydWarshallGTest, testGetDistanceDirectedGraphWithNegativeSelfLoop) {
    auto graph = this->directedGraphWithNegativeSelfLoop();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    constexpr edgeweight expectedDistance{-std::numeric_limits<edgeweight>::infinity()};
    EXPECT_THAT(this->distancesFromGetDistance(test, graph.numberOfNodes()),
                testing::Each(testing::Each(expectedDistance)));
}

TYPED_TEST_P(FloydWarshallGTest, testIsNodeInNegativeCycleDirectedGraphWithNegativeSelfLoop) {
    auto graph = this->directedGraphWithNegativeSelfLoop();
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    EXPECT_THAT(this->negativeCycleFlags(test, graph.numberOfNodes()), testing::Each(true));
}

TYPED_TEST_P(FloydWarshallGTest, testMultipleShortestDistancePaths) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    GraphT graph(11, true);
    graph.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{1});
    graph.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{1});
    graph.addEdge(NodeT{2}, NodeT{3}, EdgeWeightT{1});
    graph.addEdge(NodeT{3}, NodeT{10}, EdgeWeightT{2});
    graph.addEdge(NodeT{0}, NodeT{4}, EdgeWeightT{1});
    graph.addEdge(NodeT{4}, NodeT{5}, EdgeWeightT{1});
    graph.addEdge(NodeT{5}, NodeT{10}, EdgeWeightT{3});
    graph.addEdge(NodeT{0}, NodeT{6}, EdgeWeightT{1});
    graph.addEdge(NodeT{6}, NodeT{7}, EdgeWeightT{1});
    graph.addEdge(NodeT{7}, NodeT{8}, EdgeWeightT{1});
    graph.addEdge(NodeT{8}, NodeT{9}, EdgeWeightT{1});
    graph.addEdge(NodeT{9}, NodeT{10}, EdgeWeightT{1});

    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    constexpr edgeweight expectedShortestDistance = 5.0;
    const std::vector<NodeT> expectedPath{NodeT{0}, NodeT{4}, NodeT{5}, NodeT{10}};
    EXPECT_EQ(test.getDistance(NodeT{0}, NodeT{10}), expectedShortestDistance);
    EXPECT_THAT(test.getNodesOnShortestPath(NodeT{0}, NodeT{10}),
                testing::ElementsAreArray(expectedPath));
}

TYPED_TEST_P(FloydWarshallGTest, testGetDistances) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    GraphT graph(3, true);
    graph.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{2});
    graph.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{3});
    typename TestFixture::FloydWarshallT test(graph);
    test.run();
    EXPECT_THAT(test.getDistances(), testing::ElementsAre(testing::ElementsAre(0.0, 2.0, 5.0),
                                                          testing::ElementsAre(2.0, 0.0, 3.0),
                                                          testing::ElementsAre(5.0, 3.0, 0.0)));
}

TYPED_TEST_P(FloydWarshallGTest, testGetDistancesThrows) {
    using GraphT = typename TestFixture::GraphT;
    using FloydWarshallT = typename TestFixture::FloydWarshallT;

    GraphT graph(1, true);
    FloydWarshallT test(graph);
    try {
        test.getDistances();
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Error, run must be called first");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

REGISTER_TYPED_TEST_SUITE_P(
    FloydWarshallGTest, testConstructorThrowsUnweightedGraph, testGetDistanceThrows,
    testIsNodeInNegativeCycleThrows, testGetNodesOnShortestPathThrows,
    testGetDistanceCompleteGraphK3, testIsNodeInNegativeCycleCompleteGraphK3,
    getNodesOnShortestPathCompleteGraphK3, testGetDistanceUndirectedGraphWithNegativeEdge,
    testIsNodeInNegativeCycleUndirectedGraphWithNegativeEdge,
    getNodesOnShortestPathUndirectedGraphWithNegativeEdge,
    testGetDistanceCompleteGraphK3NegativeEdge, testIsNodeInNegativeEdgeCompleteGraphK3NegativeEdge,
    testGetNodesOnShortestPathCompleteGraphK3NegativeEdge, testGetDistanceDisconnectedGraph,
    testGetNodesOnShortestPathDisconnectedGraph, testIsNodeInNegativeCycleDisconnectedGraph,
    testGetDistanceDirectedGraphWithNegativeSelfLoop,
    testIsNodeInNegativeCycleDirectedGraphWithNegativeSelfLoop, testMultipleShortestDistancePaths,
    testGetDistances, testGetDistancesThrows);

using FloydWarshallTestTypes =
    ::testing::Types<FloydWarshallConfig<node, edgeweight>, FloydWarshallConfig<int, float>,
                     FloydWarshallConfig<int, int>>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestFloydWarshall, FloydWarshallGTest, FloydWarshallTestTypes, );

} // namespace
} // namespace NetworKit
