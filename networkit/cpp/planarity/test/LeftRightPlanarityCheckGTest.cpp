/*
 * GraphGTest.cpp
 *
 *  Created on: 05.01.2025
 *      Author: Andreas Scharf (andreas.b.scharf@gmail.com)
 */
#include <stdexcept>

#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/planarity/LeftRightPlanarityCheck.hpp>

namespace NetworKit {
namespace {

template <class NodeT_, class EdgeWeightT_, bool Weighted>
struct PlanarityGraphConfig {
    using NodeT = NodeT_;
    using EdgeWeightT = EdgeWeightT_;
    static constexpr bool weighted = Weighted;
};

template <class TestT>
class GenericLeftRightPlanarityCheckGTest : public testing::Test {
public:
    using NodeT = typename TestT::NodeT;
    using EdgeWeightT = typename TestT::EdgeWeightT;
    using GraphType = AdjListGraph<NodeT, EdgeWeightT>;
    using PlanarityCheck = GenericLeftRightPlanarityCheck<GraphType>;

    static constexpr int maxNumberOfNodes{10};

    GraphType graphWithNodes(NodeT numNodes) {
        return GraphType(numNodes, TestT::weighted, false, true);
    }

    GraphType pathGraph(NodeT numNodes) {
        GraphType graph = graphWithNodes(numNodes);
        for (NodeT i = 0; i < numNodes - 1; ++i) {
            graph.addEdge(i, i + 1);
        }
        return graph;
    }

    GraphType cycleGraph(NodeT numNodes) {
        GraphType graph = graphWithNodes(numNodes);
        for (NodeT i = 0; i < numNodes - 1; ++i) {
            graph.addEdge(i, i + 1);
        }
        if (numNodes > 2)
            graph.addEdge(numNodes - 2, NodeT{0});
        return graph;
    }

    GraphType starGraph(NodeT numNodes) {
        GraphType graph = graphWithNodes(numNodes);
        for (NodeT i = 0; i < numNodes - 1; ++i) {
            graph.addEdge(i, i + 1);
        }
        if (numNodes > 2)
            graph.addEdge(numNodes - 2, NodeT{0});
        return graph;
    }

    GraphType binaryTreeGraph(NodeT numNodes) {
        GraphType graph = graphWithNodes(numNodes);
        for (NodeT i = 0; i < numNodes; ++i) {
            const NodeT leftChild = 2 * i + 1;
            const NodeT rightChild = 2 * i + 2;
            if (leftChild < numNodes) {
                graph.addEdge(i, leftChild, static_cast<EdgeWeightT>(i));
            }
            if (rightChild < numNodes) {
                graph.addEdge(i, rightChild, static_cast<EdgeWeightT>(i));
            }
        }
        return graph;
    }

    GraphType wheelGraph(NodeT numNodes) {
        GraphType graph = graphWithNodes(numNodes);
        if (numNodes < 4) {
            throw std::invalid_argument("A wheel graph requires at least 4 nodes.");
        }
        // Form cycle
        for (NodeT i = 1; i < numNodes - 1; ++i) {
            graph.addEdge(i, i + 1);
        }
        graph.addEdge(numNodes - 1, NodeT{1}); // Close the cycle

        // Connect center to cycle
        for (NodeT i = 1; i < numNodes; ++i) {
            graph.addEdge(NodeT{0}, i);
        }
        return graph;
    }

    GraphType completeGraph(NodeT numNodes) {
        GraphType graph = graphWithNodes(numNodes);

        for (NodeT i = 0; i < numNodes; ++i) {
            for (NodeT j = i + 1; j < numNodes; ++j) {
                graph.addEdge(i, j, static_cast<EdgeWeightT>(j * (j + 1)));
            }
        }
        return graph;
    }

    GraphType gridGraph(NodeT rows, NodeT columns) {
        GraphType graph = graphWithNodes(rows * columns);
        for (NodeT row = 0; row < rows; ++row) {
            for (NodeT col = 0; col < columns; ++col) {
                const NodeT currentNode = row * columns + col;

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

    GraphType petersenGraph(NodeT n, NodeT k) {
        GraphType graph = graphWithNodes(2 * n);

        for (NodeT i = 0; i < n; ++i) {
            graph.addEdge(i, (i + 1) % n);
        }

        for (NodeT i = 0; i < n; ++i) {
            graph.addEdge(n + i, n + (i + k) % n);
        }

        for (NodeT i = 0; i < n; ++i) {
            graph.addEdge(i, n + i);
        }
        return graph;
    }
};

using PlanarityTestTypes = ::testing::Types<
    PlanarityGraphConfig<node, edgeweight, false>, PlanarityGraphConfig<node, edgeweight, true>,
    PlanarityGraphConfig<int, float, false>, PlanarityGraphConfig<int, float, true>,
    PlanarityGraphConfig<int, int, false>, PlanarityGraphConfig<int, int, true>>;

TYPED_TEST_SUITE(GenericLeftRightPlanarityCheckGTest, PlanarityTestTypes, );

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testIsPlanarThrowsIfRunIsNotCalled) {
    typename TestFixture::GraphType graph = this->graphWithNodes(0);
    typename TestFixture::PlanarityCheck test(graph);
    try {
        test.isPlanar();
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Error, run must be called first");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testNoEdgesIndexedGraphThrows) {
    typename TestFixture::GraphType graph(0, TypeParam::weighted, false, false);
    try {
        typename TestFixture::PlanarityCheck test(graph);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "The graph has no edge IDs.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testDirectedGraphThrows) {
    typename TestFixture::GraphType graph(0, TypeParam::weighted, true, true);
    try {
        typename TestFixture::PlanarityCheck test(graph);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "The graph is not an undirected graph.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testPlanarEmptyGraph) {
    typename TestFixture::GraphType graph = this->graphWithNodes(0);
    typename TestFixture::PlanarityCheck test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testPlanarSingleNode) {
    typename TestFixture::GraphType graph = this->graphWithNodes(1);
    typename TestFixture::PlanarityCheck test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testPlanarPathGraphs) {
    using NodeT = typename TestFixture::NodeT;

    for (NodeT numberOfNodes = 2; numberOfNodes <= TestFixture::maxNumberOfNodes; ++numberOfNodes) {
        typename TestFixture::GraphType graph = this->pathGraph(numberOfNodes);
        typename TestFixture::PlanarityCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testPlanarCycleGraphs) {
    using NodeT = typename TestFixture::NodeT;

    for (NodeT numberOfNodes = 2; numberOfNodes <= TestFixture::maxNumberOfNodes; ++numberOfNodes) {
        typename TestFixture::GraphType graph = this->cycleGraph(numberOfNodes);
        typename TestFixture::PlanarityCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testPlanarStarGraphs) {
    using NodeT = typename TestFixture::NodeT;

    for (NodeT numberOfNodes = 2; numberOfNodes <= TestFixture::maxNumberOfNodes; ++numberOfNodes) {
        typename TestFixture::GraphType graph = this->starGraph(numberOfNodes);
        typename TestFixture::PlanarityCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testPlanarTreeGraphs) {
    using NodeT = typename TestFixture::NodeT;

    for (NodeT numberOfNodes = 2; numberOfNodes <= TestFixture::maxNumberOfNodes; ++numberOfNodes) {
        typename TestFixture::GraphType graph = this->binaryTreeGraph(numberOfNodes);
        typename TestFixture::PlanarityCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testPlanarWheelGraphs) {
    using NodeT = typename TestFixture::NodeT;

    for (NodeT numberOfNodes = 4; numberOfNodes <= TestFixture::maxNumberOfNodes; ++numberOfNodes) {
        typename TestFixture::GraphType graph = this->wheelGraph(numberOfNodes);
        typename TestFixture::PlanarityCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testPlanarCompleteGraphs) {
    using NodeT = typename TestFixture::NodeT;

    constexpr NodeT maxNumberPlanar{5};
    for (NodeT numberOfNodes = 2; numberOfNodes < maxNumberPlanar; ++numberOfNodes) {
        typename TestFixture::GraphType graph = this->completeGraph(numberOfNodes);
        typename TestFixture::PlanarityCheck test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testNonPlanarCompleteGraphsEulerCriterium) {
    using NodeT = typename TestFixture::NodeT;

    for (NodeT numberOfNodes = 5; numberOfNodes <= TestFixture::maxNumberOfNodes; ++numberOfNodes) {
        typename TestFixture::GraphType graph = this->completeGraph(numberOfNodes);
        typename TestFixture::PlanarityCheck test(graph);
        test.run();
        EXPECT_FALSE(test.isPlanar());
    }
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testPlanarGridGraphs) {
    using NodeT = typename TestFixture::NodeT;

    for (NodeT numberOfRows = 2; numberOfRows < TestFixture::maxNumberOfNodes / 2; ++numberOfRows) {
        for (NodeT numberOfColumns = 2; numberOfColumns < TestFixture::maxNumberOfNodes / 2;
             ++numberOfColumns) {
            typename TestFixture::GraphType graph = this->gridGraph(numberOfRows, numberOfColumns);
            typename TestFixture::PlanarityCheck test(graph);
            test.run();
            EXPECT_TRUE(test.isPlanar());
        }
    }
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testNonPlanarCompleteBipartiteGraphK3_3) {
    using NodeT = typename TestFixture::NodeT;

    typename TestFixture::GraphType graph = this->graphWithNodes(6);
    graph.addEdge(NodeT{0}, NodeT{3});
    graph.addEdge(NodeT{0}, NodeT{4});
    graph.addEdge(NodeT{0}, NodeT{5});
    graph.addEdge(NodeT{1}, NodeT{3});
    graph.addEdge(NodeT{1}, NodeT{4});
    graph.addEdge(NodeT{1}, NodeT{5});
    graph.addEdge(NodeT{2}, NodeT{3});
    graph.addEdge(NodeT{2}, NodeT{4});
    graph.addEdge(NodeT{2}, NodeT{5});

    typename TestFixture::PlanarityCheck test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testNonPlanarCompleteTripartiteGraphK3_3_3) {
    using NodeT = typename TestFixture::NodeT;

    typename TestFixture::GraphType graph = this->graphWithNodes(9);
    graph.addEdge(NodeT{0}, NodeT{3});
    graph.addEdge(NodeT{0}, NodeT{4});
    graph.addEdge(NodeT{0}, NodeT{5});
    graph.addEdge(NodeT{1}, NodeT{3});
    graph.addEdge(NodeT{1}, NodeT{4});
    graph.addEdge(NodeT{1}, NodeT{5});
    graph.addEdge(NodeT{2}, NodeT{3});
    graph.addEdge(NodeT{2}, NodeT{4});
    graph.addEdge(NodeT{2}, NodeT{5});

    graph.addEdge(NodeT{0}, NodeT{6});
    graph.addEdge(NodeT{0}, NodeT{7});
    graph.addEdge(NodeT{0}, NodeT{8});
    graph.addEdge(NodeT{1}, NodeT{6});
    graph.addEdge(NodeT{1}, NodeT{7});
    graph.addEdge(NodeT{1}, NodeT{8});
    graph.addEdge(NodeT{2}, NodeT{6});
    graph.addEdge(NodeT{2}, NodeT{7});
    graph.addEdge(NodeT{2}, NodeT{8});

    graph.addEdge(NodeT{3}, NodeT{6});
    graph.addEdge(NodeT{3}, NodeT{7});
    graph.addEdge(NodeT{3}, NodeT{8});
    graph.addEdge(NodeT{4}, NodeT{6});
    graph.addEdge(NodeT{4}, NodeT{7});
    graph.addEdge(NodeT{4}, NodeT{8});
    graph.addEdge(NodeT{5}, NodeT{6});
    graph.addEdge(NodeT{5}, NodeT{7});
    graph.addEdge(NodeT{5}, NodeT{8});

    typename TestFixture::PlanarityCheck test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testOnePlanarOneNonPlanarSubGraph) {
    using NodeT = typename TestFixture::NodeT;

    typename TestFixture::GraphType graph = this->graphWithNodes(10);
    // complete bipartite graph K3,3 (non-planar)
    graph.addEdge(NodeT{0}, NodeT{3});
    graph.addEdge(NodeT{0}, NodeT{4});
    graph.addEdge(NodeT{0}, NodeT{5});
    graph.addEdge(NodeT{1}, NodeT{3});
    graph.addEdge(NodeT{1}, NodeT{4});
    graph.addEdge(NodeT{1}, NodeT{5});
    graph.addEdge(NodeT{2}, NodeT{3});
    graph.addEdge(NodeT{2}, NodeT{4});
    graph.addEdge(NodeT{2}, NodeT{5});
    // Simple cycle graph (planar)
    graph.addEdge(NodeT{6}, NodeT{7});
    graph.addEdge(NodeT{7}, NodeT{8});
    graph.addEdge(NodeT{8}, NodeT{9});
    graph.addEdge(NodeT{9}, NodeT{6});

    typename TestFixture::PlanarityCheck test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testPlanarPetersenGraphs) {
    using NodeT = typename TestFixture::NodeT;

    for (NodeT n = 3; n < TestFixture::maxNumberOfNodes; ++n) {
        for (NodeT k = 1; k <= n / 2; ++k) {
            const bool isPlanarPetersenGraph = k == 1 || (k == 2 && !(n & 1));
            if (isPlanarPetersenGraph) {
                typename TestFixture::GraphType graph = this->petersenGraph(n, k);
                typename TestFixture::PlanarityCheck test(graph);
                test.run();
                EXPECT_TRUE(test.isPlanar());
            }
        }
    }
}

TYPED_TEST(GenericLeftRightPlanarityCheckGTest, testNonPlanarPetersenGraphs) {
    using NodeT = typename TestFixture::NodeT;

    for (NodeT n = 3; n < TestFixture::maxNumberOfNodes; ++n) {
        for (NodeT k = 1; k <= n / 2; ++k) {
            const bool isNonPlanarPetersenGraph = !(k == 1 || (k == 2 && !(n & 1)));
            if (isNonPlanarPetersenGraph) {
                typename TestFixture::GraphType graph = this->petersenGraph(n, k);
                typename TestFixture::PlanarityCheck test(graph);
                test.run();
                EXPECT_FALSE(test.isPlanar());
            }
        }
    }
}

TEST(LeftRightPlanarityCheckGTest, testPlanar4eltGraph) {
    METISGraphReader reader;
    Graph graph = reader.read("input/4elt.graph");
    graph.indexEdges();
    LeftRightPlanarityCheck test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TEST(LeftRightPlanarityCheckGTest, testNonPlanarHepthGraph) {
    METISGraphReader reader;
    Graph graph = reader.read("input/hep-th.graph");
    graph.indexEdges();
    LeftRightPlanarityCheck test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

TEST(LeftRightPlanarityCheckGTest, testPlanarAirfoil1Graph) {
    METISGraphReader reader;
    Graph graph = reader.read("input/airfoil1.graph");
    graph.indexEdges();
    LeftRightPlanarityCheck test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TEST(LeftRightPlanarityCheckGTest, testNonPlanarAstroPhGraph) {
    METISGraphReader reader;
    Graph graph = reader.read("input/astro-ph.graph");
    graph.indexEdges();
    LeftRightPlanarityCheck test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

} // namespace
} // namespace NetworKit
