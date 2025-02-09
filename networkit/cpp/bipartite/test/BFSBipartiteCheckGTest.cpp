/*
 * BFSBipartiteCheckGTest.cpp
 *
 *  Created on: 09.02.2025
 *      Author: Andreas Scharf (andreas.b.scharf@gmail.com)
 */

#include <gtest/gtest.h>
#include <networkit/bipartite/BFSBipartiteCheck.hpp>
#include <networkit/io/METISGraphReader.hpp>
namespace NetworKit {

class BFSBipartiteCheckGTest : public testing::Test {
public:
    static constexpr int maxNumberOfNodes{10};
    Graph binaryTreeGraph(count numNodes) {
        Graph graph(numNodes, true, false, true);
        for (count i = 0; i < numNodes; ++i) {
            count leftChild = 2 * i + 1;
            count rightChild = 2 * i + 2;
            if (leftChild < numNodes) {
                graph.addEdge(i, leftChild, static_cast<double>(i));
            }
            if (rightChild < numNodes) {
                graph.addEdge(i, rightChild, static_cast<double>(i));
            }
        }
        return graph;
    }

    Graph completeGraph(count numNodes) {
        Graph graph(numNodes, true);

        for (count i = 0; i < numNodes; ++i) {
            for (count j = i + 1; j < numNodes; ++j) {
                graph.addEdge(i, j, static_cast<double>(j * (j + 1)));
            }
        }
        return graph;
    }
};

TEST_F(BFSBipartiteCheckGTest, testDirectedGraphThrows) {
    Graph graph(0, false, true, false);
    try {
        BFSBipartiteCheck test(graph);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "The graph is not an undirected graph.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(BFSBipartiteCheckGTest, testIsBipartiteThrowsIfRunIsNotCalled) {
    Graph graph{};
    BFSBipartiteCheck test(graph);
    try {
        test.isBipartite();
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Error, run must be called first");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(BFSBipartiteCheckGTest, testBipartiteEmptyGraph) {
    Graph graph{};
    BFSBipartiteCheck test(graph);
    test.run();
    EXPECT_TRUE(test.isBipartite());
}

TEST_F(BFSBipartiteCheckGTest, testBipartiteSingleNodesGraphs) {
    for (count i = 0; i < maxNumberOfNodes; ++i) {
        Graph graph{};
        graph.addNodes(1);
        BFSBipartiteCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isBipartite());
    }
}

TEST_F(BFSBipartiteCheckGTest, testBipartiteBinaryTreeGraphs) {
    for (count numberOfNodes = 2; numberOfNodes <= maxNumberOfNodes; ++numberOfNodes) {
        Graph graph = binaryTreeGraph(numberOfNodes);
        BFSBipartiteCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isBipartite());
    }
}

TEST_F(BFSBipartiteCheckGTest, testNonBipartiteCompleteGraphs) {
    for (count numberOfNodes = 3; numberOfNodes <= maxNumberOfNodes; ++numberOfNodes) {
        Graph graph = completeGraph(numberOfNodes);
        BFSBipartiteCheck test(graph);
        test.run();
        EXPECT_FALSE(test.isBipartite());
    }
}

TEST_F(BFSBipartiteCheckGTest, testBipartiteGrid5x5DistArchGraph) {
    METISGraphReader reader;
    Graph graph = reader.read("input/grid-5x5-dist-arch.graph");
    BFSBipartiteCheck test(graph);
    test.run();
    EXPECT_TRUE(test.isBipartite());
}

TEST_F(BFSBipartiteCheckGTest, testNonBipartiteAirfoil1Graph) {
    METISGraphReader reader;
    Graph graph = reader.read("input/airfoil1.graph");
    BFSBipartiteCheck test(graph);
    test.run();
    EXPECT_FALSE(test.isBipartite());
}

TEST_F(BFSBipartiteCheckGTest, testNonBipartiteHepThGraph) {
    METISGraphReader reader;
    Graph graph = reader.read("input/hep-th.graph");
    BFSBipartiteCheck test(graph);
    test.run();
    EXPECT_FALSE(test.isBipartite());
}

TEST_F(BFSBipartiteCheckGTest, testCompleteBipartiteGraphK3_3) {
    Graph graph(6);
    graph.addEdge(0, 3);
    graph.addEdge(0, 4);
    graph.addEdge(0, 5);
    graph.addEdge(1, 3);
    graph.addEdge(1, 4);
    graph.addEdge(1, 5);
    graph.addEdge(2, 3);
    graph.addEdge(2, 4);
    graph.addEdge(2, 5);

    BFSBipartiteCheck test(graph);
    test.run();
    EXPECT_TRUE(test.isBipartite());
}

} // namespace NetworKit
