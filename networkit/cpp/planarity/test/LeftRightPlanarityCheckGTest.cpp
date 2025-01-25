/*
 * GraphGTest.cpp
 *
 *  Created on: 05.01.2025
 *      Author: Andreas Scharf (andreas.b.scharf@gmail.com)
 */
#include <gtest/gtest.h>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/planarity/LeftRightPlanarityCheck.hpp>

namespace NetworKit {
class LeftRightPlanarityCheckGTest : public testing::Test {
public:
    LeftRightPlanarityCheckGTest() = default;
    static constexpr int maxNumberOfNodes{10};
    Graph pathGraph(count numNodes) {
        Graph graph;
        graph.addNodes(numNodes);
        for (count i = 0; i < numNodes - 1; ++i) {
            graph.addEdge(i, i + 1);
        }
        return graph;
    }

    Graph cycleGraph(count numNodes) {
        Graph graph;
        graph.addNodes(numNodes);
        for (count i = 0; i < numNodes - 1; ++i) {
            graph.addEdge(i, i + 1);
        }
        if (numNodes > 2)
            graph.addEdge(numNodes - 2, 0);
        return graph;
    }

    Graph starGraph(count numNodes) {
        Graph graph;
        graph.addNodes(numNodes);
        for (count i = 0; i < numNodes - 1; ++i) {
            graph.addEdge(i, i + 1);
        }
        if (numNodes > 2)
            graph.addEdge(numNodes - 2, 0);
        return graph;
    }

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

    Graph wheelGraph(count numNodes) {
        Graph graph(numNodes, false, false, true);
        if (numNodes < 4) {
            throw std::invalid_argument("A wheel graph requires at least 4 nodes.");
        }
        // Form cycle
        for (count i = 1; i < numNodes - 1; ++i) {
            graph.addEdge(i, i + 1);
        }
        graph.addEdge(numNodes - 1, 1); // Close the cycle

        // Connect center to cycle
        for (count i = 1; i < numNodes; ++i) {
            graph.addEdge(0, i);
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

    Graph gridGraph(count rows, count columns) {
        Graph graph(rows * columns);
        for (count row = 0; row < rows; ++row) {
            for (count col = 0; col < columns; ++col) {
                count currentNode = row * columns + col;

                // Connect to the right neighbor
                if (col + 1 < columns) {
                    graph.addEdge(currentNode, currentNode + 1);
                }

                // Connect to the bottom neighbor
                if (row + 1 < rows) {
                    graph.addEdge(currentNode, currentNode + columns);
                }
            }
        }
        return graph;
    }

    Graph petersenGraph(count n, count k) {
        Graph graph(2 * n);

        for (count i = 0; i < n; ++i) {
            graph.addEdge(i, (i + 1) % n);
        }

        for (count i = 0; i < n; ++i) {
            graph.addEdge(n + i, n + (i + k) % n);
        }

        for (count i = 0; i < n; ++i) {
            graph.addEdge(i, n + i);
        }

        return graph;
    }
};

TEST_F(LeftRightPlanarityCheckGTest, testDirectedGraphThrows) {
    Graph graph(0, false, true, false);
    try {
        LeftRightPlanarityCheck test(graph);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "The graph is not an undirected graph.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(LeftRightPlanarityCheckGTest, testIsPlanarThrowsIfRunIsNotCalled) {
    Graph graph{};
    LeftRightPlanarityCheck test(graph);
    try {
        test.isPlanar();
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Error, run must be called first");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(LeftRightPlanarityCheckGTest, testPlanarEmptyGraph) {
    Graph graph{};
    LeftRightPlanarityCheck test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TEST_F(LeftRightPlanarityCheckGTest, testPlanarSingleNode) {
    Graph graph{1};
    LeftRightPlanarityCheck test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TEST_F(LeftRightPlanarityCheckGTest, testPlanarPathGraphs) {
    for (count numberOfNodes = 2; numberOfNodes <= maxNumberOfNodes; ++numberOfNodes) {
        Graph graph = pathGraph(numberOfNodes);
        LeftRightPlanarityCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityCheckGTest, testPlanarCycleGraphs) {
    for (count numberOfNodes = 2; numberOfNodes <= maxNumberOfNodes; ++numberOfNodes) {
        Graph graph = cycleGraph(numberOfNodes);
        LeftRightPlanarityCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityCheckGTest, testPlanarStarGraphs) {
    for (count numberOfNodes = 2; numberOfNodes <= maxNumberOfNodes; ++numberOfNodes) {
        Graph graph = starGraph(numberOfNodes);
        LeftRightPlanarityCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityCheckGTest, testPlanarTreeGraphs) {
    for (count numberOfNodes = 2; numberOfNodes <= maxNumberOfNodes; ++numberOfNodes) {
        Graph graph = binaryTreeGraph(numberOfNodes);
        LeftRightPlanarityCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityCheckGTest, testPlanarWheelGraphs) {
    for (count numberOfNodes = 4; numberOfNodes <= maxNumberOfNodes; ++numberOfNodes) {
        Graph graph = wheelGraph(numberOfNodes);
        LeftRightPlanarityCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityCheckGTest, testPlanarCompleteGraphs) {
    constexpr count maxNumberPlanar{5};
    for (count numberOfNodes = 2; numberOfNodes < maxNumberPlanar; ++numberOfNodes) {
        Graph graph = completeGraph(numberOfNodes);
        LeftRightPlanarityCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityCheckGTest, testNonPlanarCompleteGraphsEulerCriterium) {
    for (count numberOfNodes = 5; numberOfNodes <= maxNumberOfNodes; ++numberOfNodes) {
        Graph graph = completeGraph(numberOfNodes);
        LeftRightPlanarityCheck test(graph);
        test.run();
        EXPECT_FALSE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityCheckGTest, testPlanarGridGraphs) {

    for (count numberOfRows = 2; numberOfRows < maxNumberOfNodes / 2; ++numberOfRows) {
        for (count numberOfColumns = 2; numberOfColumns < maxNumberOfNodes / 2; ++numberOfColumns) {
            Graph graph = gridGraph(numberOfRows, numberOfColumns);
            LeftRightPlanarityCheck test(graph);
            test.run();
            EXPECT_TRUE(test.isPlanar());
        }
    }
}

TEST_F(LeftRightPlanarityCheckGTest, testNonPlanarCompleteBipartiteGraphK3_3) {
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

    LeftRightPlanarityCheck test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

TEST_F(LeftRightPlanarityCheckGTest, testNonPlanarCompleteTripartiteGraphK3_3_3) {
    Graph graph(9);
    graph.addEdge(0, 3);
    graph.addEdge(0, 4);
    graph.addEdge(0, 5);
    graph.addEdge(1, 3);
    graph.addEdge(1, 4);
    graph.addEdge(1, 5);
    graph.addEdge(2, 3);
    graph.addEdge(2, 4);
    graph.addEdge(2, 5);

    graph.addEdge(0, 6);
    graph.addEdge(0, 7);
    graph.addEdge(0, 8);
    graph.addEdge(1, 6);
    graph.addEdge(1, 7);
    graph.addEdge(1, 8);
    graph.addEdge(2, 6);
    graph.addEdge(2, 7);
    graph.addEdge(2, 8);

    graph.addEdge(3, 6);
    graph.addEdge(3, 7);
    graph.addEdge(3, 8);
    graph.addEdge(4, 6);
    graph.addEdge(4, 7);
    graph.addEdge(4, 8);
    graph.addEdge(5, 6);
    graph.addEdge(5, 7);
    graph.addEdge(5, 8);
    LeftRightPlanarityCheck test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

TEST_F(LeftRightPlanarityCheckGTest, testOnePlanarOneNonPlanarSubGraph) {
    Graph graph(10);
    // complete bipartite graph K3,3 (non-planar)
    graph.addEdge(0, 3);
    graph.addEdge(0, 4);
    graph.addEdge(0, 5);
    graph.addEdge(1, 3);
    graph.addEdge(1, 4);
    graph.addEdge(1, 5);
    graph.addEdge(2, 3);
    graph.addEdge(2, 4);
    graph.addEdge(2, 5);
    // Simple cycle graph (planar)
    graph.addEdge(6, 7);
    graph.addEdge(7, 8);
    graph.addEdge(8, 9);
    graph.addEdge(9, 6);
    LeftRightPlanarityCheck test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

TEST_F(LeftRightPlanarityCheckGTest, testPlanarPetersenGraphs) {
    for (count n = 3; n < maxNumberOfNodes; ++n) {
        for (count k = 1; k <= std::floor(n / 2); ++k) {
            const bool isPlanarPetersenGraph = k == 1 || (k == 2 && !(n & 1));
            if (isPlanarPetersenGraph) {
                Graph graph = petersenGraph(n, k);
                LeftRightPlanarityCheck test(graph);
                test.run();
                EXPECT_TRUE(test.isPlanar());
            }
        }
    }
}

TEST_F(LeftRightPlanarityCheckGTest, testNonPlanarPetersenGraphs) {
    for (count n = 3; n < maxNumberOfNodes; ++n) {
        for (count k = 1; k <= std::floor(n / 2); ++k) {
            const bool isNonPlanarPetersenGraph = !(k == 1 || (k == 2 && !(n & 1)));
            if (isNonPlanarPetersenGraph) {
                Graph graph = petersenGraph(n, k);
                LeftRightPlanarityCheck test(graph);
                test.run();
                EXPECT_FALSE(test.isPlanar());
            }
        }
    }
}

TEST_F(LeftRightPlanarityCheckGTest, testPlanar4eltGraph) {
    METISGraphReader reader;
    Graph graph = reader.read("input/4elt.graph");
    LeftRightPlanarityCheck test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TEST_F(LeftRightPlanarityCheckGTest, testNonPlanarHepthGraph) {
    METISGraphReader reader;
    Graph graph = reader.read("input/hep-th.graph");
    LeftRightPlanarityCheck test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

TEST_F(LeftRightPlanarityCheckGTest, testPlanarAirfoil1Graph) {
    METISGraphReader reader;
    Graph graph = reader.read("input/airfoil1.graph");
    LeftRightPlanarityCheck test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TEST_F(LeftRightPlanarityCheckGTest, testNonPlanarAstroPhGraph) {
    METISGraphReader reader;
    Graph graph = reader.read("input/astro-ph.graph");
    LeftRightPlanarityCheck test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

} // namespace NetworKit
