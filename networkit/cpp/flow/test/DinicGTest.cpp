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

TEST_F(DinicGTest, ThreeDisjointPaths) {
    Graph G(5, true, true);
    G.addEdge(0,1,1); G.addEdge(1,4,1);
    G.addEdge(0,2,1); G.addEdge(2,4,1);
    G.addEdge(0,3,1); G.addEdge(3,4,1);

    Dinic test(G, 0, 4);
    test.run();
    EXPECT_DOUBLE_EQ(3.0, test.getMaxFlow());
}

TEST_F(DinicGTest, testThreeCycleWithTailGraph) {
    Graph graph(4, true, true);
    graph.addEdge(0, 1, 0.3);
    graph.addEdge(1, 2, 0.6);
    graph.addEdge(2, 0, 0.9);
    graph.addEdge(2, 3, 0.7);
    constexpr double expectedFlow0to3{0.3};
    Dinic test0to3(graph, 0, 3);
    test0to3.run();
    EXPECT_DOUBLE_EQ(expectedFlow0to3, test0to3.getMaxFlow());

    constexpr double expectedFlow1to3{0.6};
    Dinic test1to3(graph, 1, 3);
    test1to3.run();
    EXPECT_DOUBLE_EQ(expectedFlow1to3, test1to3.getMaxFlow());
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
    constexpr double expectedFlow0to7{3.0};
    Dinic test0to7(graph, 0, 7);
    test0to7.run();
    EXPECT_DOUBLE_EQ(expectedFlow0to7, test0to7.getMaxFlow());

    constexpr double expectedFlow3to7{2.0};
    Dinic test3to7(graph, 3, 7);
    test3to7.run();
    EXPECT_DOUBLE_EQ(expectedFlow3to7, test3to7.getMaxFlow());

    constexpr double expectedFlow0to5{2.0};
    Dinic test0to5(graph, 0, 5);
    test0to5.run();
    EXPECT_DOUBLE_EQ(expectedFlow0to5, test0to5.getMaxFlow());

    constexpr double expectedFlow2to4{1.0};
    Dinic test2to4(graph, 2, 4);
    test2to4.run();
    EXPECT_DOUBLE_EQ(expectedFlow2to4, test2to4.getMaxFlow());
}

TEST_F(DinicGTest, testMaxFlowDiamondWithCrossGraph) {
    Graph graph(4, true, true);
    graph.addEdge(0, 1, 10.0);
    graph.addEdge(0, 2, 10.0);
    graph.addEdge(1, 2, 5.0);
    graph.addEdge(1, 3, 10.0);
    graph.addEdge(2, 3, 10.0);

    constexpr double expectedFlow0to3{20.0};
    Dinic test0to3(graph, 0, 3);
    test0to3.run();
    EXPECT_DOUBLE_EQ(expectedFlow0to3, test0to3.getMaxFlow());

    constexpr double expectedFlow0to2{15.0};
    Dinic test0to2(graph, 0, 2);
    test0to2.run();
    EXPECT_DOUBLE_EQ(expectedFlow0to2, test0to2.getMaxFlow());
}

TEST_F(DinicGTest, testMaxFlowDisconnectedGraphs) {
    Graph graph(7, true, true);
    graph.addEdge(0, 1, 10.0);
    graph.addEdge(1, 2, 5.0);
    graph.addEdge(2, 3, 7.0);

    graph.addEdge(4, 5, 11.0);
    graph.addEdge(5, 6, 10.0);
    constexpr double expectedFlow{0.0};
    Dinic test(graph, 0, 5);
    test.run();
    EXPECT_DOUBLE_EQ(expectedFlow, test.getMaxFlow());
}

} // namespace NetworKit
