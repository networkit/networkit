/*  DinicGTest.cpp
 *
 *	Created on: 20.06.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <networkit/flow/Dinic.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
class DinicGTest : public testing::Test {
public:
    DinicGTest() = default;
};

TEST_F(DinicGTest, testConstructorThrowsForUndirectedGraph) {
    Graph graph(2, true, false);
    try {
        Dinic test(graph, /* source */ 0, /*target */ 1);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Dinic algorithm requires directed graph!");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(DinicGTest, testConstructorThrowsForUnweightedGraph) {
    Graph graph(2, false, true);
    try {
        Dinic test(graph, /* source */ 0, /*target */ 1);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Dinic algorithm requires weighted graph!");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(DinicGTest, testConstructorThrowsForSameSourceAndTarget) {
    Graph graph(2, true, true);
    try {
        Dinic test(graph, /* source */ 0, /*target */ 0);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(),
                     "Dinic algorithm requires `source` and `target` node to be different!");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(DinicGTest, testRunNotCalledThrows) {
    Graph graph(2, true, true);
    Dinic test(graph, /* source */ 0, /*target */ 1);
    try {
        test.getMaxFlow();
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Error, run must be called first");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(DinicGTest, testNegativeWeightsThrows) {
    Graph G(3, true, true);
    G.addEdge(0, 1, 1);
    G.addEdge(1, 2, -1);

    Dinic test(G, 0, 2);
    EXPECT_THROW(test.run(), std::runtime_error);
}

TEST_F(DinicGTest, ThreeDisjointPaths) {
    Graph G(5, true, true);
    G.addEdge(0, 1, 1);
    G.addEdge(1, 4, 1);
    G.addEdge(0, 2, 1);
    G.addEdge(2, 4, 1);
    G.addEdge(0, 3, 1);
    G.addEdge(3, 4, 1);

    Dinic test(G, 0, 4);
    test.run();
    EXPECT_DOUBLE_EQ(test.getMaxFlow(), /*expectedFlow*/ 3.0);
}

TEST_F(DinicGTest, ThreeDisjointPathsWithParallelEdges) {
    Graph G(5, true, true);
    G.addEdge(0, 1, 1);
    G.addEdge(1, 0, 1);
    G.addEdge(1, 4, 1);
    G.addEdge(0, 2, 1);
    G.addEdge(2, 0, 1);
    G.addEdge(2, 4, 1);
    G.addEdge(0, 3, 1);
    G.addEdge(3, 0, 1);
    G.addEdge(3, 4, 1);
    G.indexEdges();
    Dinic test(G, 0, 4);
    test.run();
    EXPECT_DOUBLE_EQ(test.getMaxFlow(), /*expectedFlow*/ 3.0);
}

TEST_F(DinicGTest, testThreeCycleWithTailGraph) {
    Graph graph(4, true, true);
    graph.addEdge(0, 1, 0.3);
    graph.addEdge(1, 2, 0.6);
    graph.addEdge(2, 0, 0.9);
    graph.addEdge(2, 3, 0.7);

    Dinic test0to3(graph, 0, 3);
    test0to3.run();
    EXPECT_DOUBLE_EQ(test0to3.getMaxFlow(), /*expectedFlow0to3*/ 0.3);

    Dinic test1to3(graph, 1, 3);
    test1to3.run();
    EXPECT_DOUBLE_EQ(test1to3.getMaxFlow(), /*expectedFlow1to3*/ 0.6);
}

TEST_F(DinicGTest, testThreeCycleWithTailGraphWithParallelEdges) {
    Graph graph(4, true, true);
    graph.addEdge(0, 1, 0.3);
    graph.addEdge(1, 0, 1.3);
    graph.addEdge(1, 2, 0.6);
    graph.addEdge(2, 1, 1.6);
    graph.addEdge(2, 0, 0.9);
    graph.addEdge(0, 2, 1.9);
    graph.addEdge(2, 3, 0.7);

    Dinic test0to3(graph, 0, 3);
    test0to3.run();
    EXPECT_DOUBLE_EQ(test0to3.getMaxFlow(), /*expectedFlow0to3*/ 0.7);

    Dinic test1to3(graph, 1, 3);
    test1to3.run();
    EXPECT_DOUBLE_EQ(test1to3.getMaxFlow(), /*expectedFlow1to3*/ 0.7);
}

TEST_F(DinicGTest, testFourLayeredDAG) {
    Graph graph(8, true, true);
    graph.addEdge(0, 1, 1.0);
    graph.addEdge(0, 2, 1.0);
    graph.addEdge(0, 3, 1.0);
    graph.addEdge(1, 4, 1.0);
    graph.addEdge(2, 4, 1.0);
    graph.addEdge(2, 5, 1.0);
    graph.addEdge(3, 5, 1.0);
    graph.addEdge(3, 6, 1.0);
    graph.addEdge(4, 7, 1.0);
    graph.addEdge(5, 7, 1.0);
    graph.addEdge(6, 7, 1.0);

    Dinic test0to7(graph, 0, 7);
    test0to7.run();
    EXPECT_DOUBLE_EQ(test0to7.getMaxFlow(), /*expectedFlow0to7*/ 3.0);

    Dinic test3to7(graph, 3, 7);
    test3to7.run();
    EXPECT_DOUBLE_EQ(test3to7.getMaxFlow(), /*expectedFlow3to7*/ 2.0);

    Dinic test0to5(graph, 0, 5);
    test0to5.run();
    EXPECT_DOUBLE_EQ(test0to5.getMaxFlow(), /*expectedFlow0to5*/ 2.0);

    Dinic test2to4(graph, 2, 4);
    test2to4.run();
    EXPECT_DOUBLE_EQ(test2to4.getMaxFlow(), /*expectedFlow2to4*/ 1.0);
}

TEST_F(DinicGTest, testMaxFlowDiamondWithCrossGraph) {
    Graph graph(4, true, true);
    graph.addEdge(0, 1, 10.0);
    graph.addEdge(0, 2, 10.0);
    graph.addEdge(1, 2, 5.0);
    graph.addEdge(1, 3, 10.0);
    graph.addEdge(2, 3, 10.0);

    Dinic test0to3(graph, 0, 3);
    test0to3.run();
    EXPECT_DOUBLE_EQ(test0to3.getMaxFlow(), /*expectedFlow0to3*/ 20.0);

    Dinic test0to2(graph, 0, 2);
    test0to2.run();
    EXPECT_DOUBLE_EQ(test0to2.getMaxFlow(), /*expectedFlow0to2*/ 15.0);
}

TEST_F(DinicGTest, testMaxFlowDisconnectedGraphs) {
    Graph graph(7, true, true);
    graph.addEdge(0, 1, 10.0);
    graph.addEdge(1, 2, 5.0);
    graph.addEdge(2, 3, 7.0);

    graph.addEdge(4, 5, 11.0);
    graph.addEdge(5, 6, 10.0);
    Dinic test(graph, 0, 5);
    test.run();
    EXPECT_DOUBLE_EQ(test.getMaxFlow(), /*expectedFlow*/ 0.0);
}

TEST_F(DinicGTest, testFourLayerDAGChangedIterationOrder) {
    Graph G(8, /*weighted=*/true, /*directed=*/true);
    G.addEdge(0, 2, 1.0);
    G.addEdge(0, 3, 1.0);
    G.addEdge(0, 1, 1.0);

    G.addEdge(1, 4, 1.0);
    G.addEdge(2, 4, 1.0);
    G.addEdge(2, 5, 1.0);
    G.addEdge(3, 5, 1.0);
    G.addEdge(3, 6, 1.0);

    G.addEdge(4, 7, 1.0);
    G.addEdge(5, 7, 1.0);
    G.addEdge(6, 7, 1.0);
    G.indexEdges();
    Dinic algo(G, /*src=*/0, /*dst=*/7);
    algo.run();

    EXPECT_DOUBLE_EQ(algo.getMaxFlow(), /*expectedFlow*/ 3.0);
}

TEST_F(DinicGTest, testNumericalStabilityDecimalSplits) {
    Graph G(7, /*weighted=*/true, /*directed=*/true);

    G.addEdge(0, 1, 1.0);
    G.addEdge(1, 2, 0.1);
    G.addEdge(2, 6, 0.1);
    G.addEdge(1, 3, 0.2);
    G.addEdge(3, 6, 0.2);
    G.addEdge(1, 4, 0.3);
    G.addEdge(4, 6, 0.3);
    G.addEdge(1, 5, 0.4);
    G.addEdge(5, 6, 0.4);

    // Add a sub-epsilon direct edge that should be ignored by tolerance gating.
    G.addEdge(0, 6, 1e-18);

    Dinic algo(G, 0, 6);
    algo.run();
    EXPECT_NEAR(algo.getMaxFlow(), 1.0, 1e-12);
}

TEST_F(DinicGTest, testNumericalStabilityTinyScale) {
    constexpr double scale = 1e-9;
    Graph G(7, /*weighted=*/true, /*directed=*/true);

    G.addEdge(0, 1, 1.0 * scale);
    G.addEdge(1, 2, 0.1 * scale);
    G.addEdge(2, 6, 0.1 * scale);
    G.addEdge(1, 3, 0.2 * scale);
    G.addEdge(3, 6, 0.2 * scale);
    G.addEdge(1, 4, 0.3 * scale);
    G.addEdge(4, 6, 0.3 * scale);
    G.addEdge(1, 5, 0.4 * scale);
    G.addEdge(5, 6, 0.4 * scale);

    // Add a sub-epsilon direct edge that should be ignored by tolerance gating.
    G.addEdge(0, 6, 1e-18);

    Dinic algo(G, 0, 6);
    algo.run();
    EXPECT_NEAR(algo.getMaxFlow(), 1.0 * scale, 1e-15);
}
} // namespace NetworKit
