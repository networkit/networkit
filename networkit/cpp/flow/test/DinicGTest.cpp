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

TEST_F(DinicGTest, testSimpleGraph) {
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

TEST_F(DinicGTest, testMaxFlow) {
    Graph graph(4, true, true);
    graph.addEdge(0, 1, 0.3);
    graph.addEdge(1, 2, 0.6);
    graph.addEdge(2, 0, 0.9);
    graph.addEdge(2, 3, 0.7);
    constexpr double expectedFlow{0.3};
    Dinic test(graph, 0, 3);
    test.run();
    EXPECT_DOUBLE_EQ(expectedFlow, test.getMaxFlow());
}

TEST_F(DinicGTest, testMaxFlow2) {
    Graph graph(4, true, true);
    graph.addEdge(0, 1, 0.4);
    graph.addEdge(1, 2, 0.6);
    graph.addEdge(2, 0, 0.9);
    graph.addEdge(2, 3, 0.7);
    constexpr double expectedFlow{0.4};
    Dinic test(graph, 0, 3);
    test.run();
    EXPECT_DOUBLE_EQ(expectedFlow, test.getMaxFlow());
}

} // namespace NetworKit
