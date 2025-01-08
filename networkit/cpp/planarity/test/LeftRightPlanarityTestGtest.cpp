//
// Created by andreas on 05.01.25.
//
#include <gtest/gtest.h>
#include <networkit/planarity/LeftRightPlanarityTest.hpp>

namespace NetworKit {
class LeftRightPlanarityTestGTest : public testing::Test {
public:
    LeftRightPlanarityTestGTest() = default;
    static constexpr int max_number_of_nodes{100};
    static Graph pathGraph(const count numNodes) {
        Graph graph;
        graph.addNodes(numNodes);
        for (count i{}; i < numNodes - 1; ++i) {
            graph.addEdge(i, i + 1);
        }
        return graph;
    }

    static Graph cycleGraph(const count numNodes) {
        Graph graph;
        graph.addNodes(numNodes);
        for (count i{}; i < numNodes - 1; ++i) {
            graph.addEdge(i, i + 1);
        }
        if (numNodes > 2)
            graph.addEdge(numNodes - 2, 0);
        return graph;
    }

    static Graph starGraph(const count numNodes) {
        Graph graph;
        graph.addNodes(numNodes);
        for (count i{}; i < numNodes - 1; ++i) {
            graph.addEdge(i, i + 1);
        }
        if (numNodes > 2)
            graph.addEdge(numNodes - 2, 0);
        return graph;
    }

    static Graph binaryTreeGraph(const count numNodes) {
        Graph graph(numNodes);
        for (count i{}; i < numNodes; ++i) {
            count leftChild = 2 * i + 1;
            count rightChild = 2 * i + 2;
            if (leftChild < numNodes) {
                graph.addEdge(i, leftChild);
            }
            if (rightChild < numNodes) {
                graph.addEdge(i, rightChild);
            }
        }
        return graph;
    }

    static Graph wheelGraph(const count numNodes) {
        Graph graph(numNodes);
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

    static Graph completeGraph(const count numNodes) {
        Graph graph(numNodes);

        for (count i = 0; i < numNodes; ++i) {
            for (count j = i + 1; j < numNodes; ++j) {
                graph.addEdge(i, j);
            }
        }
        return graph;
    }

    static Graph gridGraph(const count rows, const count columns) {
        Graph graph(rows * columns);
        for (count row = 0; row < rows; ++row) {
            for (count col = 0; col < columns; ++col) {
                const count currentNode = row * columns + col;

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
};

TEST_F(LeftRightPlanarityTestGTest, EmptyGraph) {

    auto graph = Graph{};
    LeftRightPlanarityTest test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TEST_F(LeftRightPlanarityTestGTest, SingleNode) {

    auto graph = Graph{1};
    LeftRightPlanarityTest test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TEST_F(LeftRightPlanarityTestGTest, PathGraph) {

    for (count numberOfNodes{2}; numberOfNodes <= max_number_of_nodes; ++numberOfNodes) {
        auto graph = pathGraph(numberOfNodes);
        LeftRightPlanarityTest test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityTestGTest, CycleGraph) {

    for (count numberOfNodes{2}; numberOfNodes <= max_number_of_nodes; ++numberOfNodes) {
        auto graph = cycleGraph(numberOfNodes);
        LeftRightPlanarityTest test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityTestGTest, StarGraph) {

    for (count numberOfNodes{2}; numberOfNodes <= max_number_of_nodes; ++numberOfNodes) {
        auto graph = starGraph(numberOfNodes);
        LeftRightPlanarityTest test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityTestGTest, TreeGraph) {

    for (count numberOfNodes{2}; numberOfNodes <= max_number_of_nodes; ++numberOfNodes) {
        auto graph = binaryTreeGraph(numberOfNodes);
        LeftRightPlanarityTest test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityTestGTest, WheelGraph) {

    for (count numberOfNodes{4}; numberOfNodes <= max_number_of_nodes; ++numberOfNodes) {
        auto graph = wheelGraph(numberOfNodes);
        LeftRightPlanarityTest test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityTestGTest, CompleteGraph) {
    for (count numberOfNodes{2}; numberOfNodes <= max_number_of_nodes; ++numberOfNodes) {
        auto graph = completeGraph(numberOfNodes);
        LeftRightPlanarityTest test(graph);
        test.run();
        if (numberOfNodes < 5)
            EXPECT_TRUE(test.isPlanar());
        else
            EXPECT_FALSE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityTestGTest, GridGraph) {

    for (count numberOfRows{2}; numberOfRows < 11; ++numberOfRows) {
        for (count numberOfColumns{2}; numberOfColumns < 11; ++numberOfColumns) {
            auto graph = gridGraph(numberOfRows, numberOfColumns);
            LeftRightPlanarityTest test(graph);
            test.run();
            EXPECT_TRUE(test.isPlanar());
        }
    }
}

TEST_F(LeftRightPlanarityTestGTest, GeneralizedPetersenGraph5_1) {
    Graph graph(10);
    graph.addEdge(0, 1);
    graph.addEdge(0, 4);
    graph.addEdge(0, 5);
    graph.addEdge(1, 2);
    graph.addEdge(1, 6);
    graph.addEdge(2, 3);
    graph.addEdge(2, 7);
    graph.addEdge(3, 4);
    graph.addEdge(3, 8);
    graph.addEdge(4, 9);
    graph.addEdge(5, 6);
    graph.addEdge(5, 9);
    graph.addEdge(6, 7);
    graph.addEdge(7, 8);
    graph.addEdge(8, 9);

    LeftRightPlanarityTest test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TEST_F(LeftRightPlanarityTestGTest, GeneralizedPetersenGraph5_2) {
    Graph graph(10);
    graph.addEdge(0, 2);
    graph.addEdge(0, 3);
    graph.addEdge(0, 5);
    graph.addEdge(1, 3);
    graph.addEdge(1, 4);
    graph.addEdge(1, 6);
    graph.addEdge(2, 4);
    graph.addEdge(2, 7);
    graph.addEdge(3, 8);
    graph.addEdge(4, 9);
    graph.addEdge(5, 6);
    graph.addEdge(5, 9);
    graph.addEdge(6, 7);
    graph.addEdge(7, 8);
    graph.addEdge(8, 9);
    LeftRightPlanarityTest test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

TEST_F(LeftRightPlanarityTestGTest, CompleteBipartiteGraphK3_3) {
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

    LeftRightPlanarityTest test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

TEST_F(LeftRightPlanarityTestGTest, CompleteTripartiteGraphK3_3_3) {
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
    LeftRightPlanarityTest test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

TEST_F(LeftRightPlanarityTestGTest, RandomGeneratedPlanarGraph50Nodes) {
    Graph graph(50);
    graph.addEdge(0, 20);
    graph.addEdge(1, 21);
    graph.addEdge(2, 14);
    graph.addEdge(2, 40);
    graph.addEdge(3, 6);
    graph.addEdge(3, 12);
    graph.addEdge(3, 15);
    graph.addEdge(4, 10);
    graph.addEdge(4, 30);
    graph.addEdge(4, 47);
    graph.addEdge(5, 7);
    graph.addEdge(5, 20);
    graph.addEdge(6, 15);
    graph.addEdge(7, 11);
    graph.addEdge(8, 10);
    graph.addEdge(8, 28);
    graph.addEdge(9, 23);
    graph.addEdge(9, 44);
    graph.addEdge(10, 11);
    graph.addEdge(10, 12);
    graph.addEdge(11, 15);
    graph.addEdge(11, 46);
    graph.addEdge(12, 16);
    graph.addEdge(12, 26);
    graph.addEdge(13, 26);
    graph.addEdge(13, 37);
    graph.addEdge(14, 19);
    graph.addEdge(15, 30);
    graph.addEdge(15, 35);
    graph.addEdge(15, 38);
    graph.addEdge(15, 39);
    graph.addEdge(16, 43);
    graph.addEdge(17, 35);
    graph.addEdge(18, 37);
    graph.addEdge(18, 43);
    graph.addEdge(20, 26);
    graph.addEdge(21, 27);
    graph.addEdge(21, 33);
    graph.addEdge(23, 25);
    graph.addEdge(24, 33);
    graph.addEdge(25, 36);
    graph.addEdge(26, 32);
    graph.addEdge(27, 42);
    graph.addEdge(28, 42);
    graph.addEdge(29, 40);
    graph.addEdge(30, 48);
    graph.addEdge(33, 34);
    graph.addEdge(34, 38);
    graph.addEdge(34, 45);
    graph.addEdge(35, 37);
    graph.addEdge(35, 46);
    graph.addEdge(36, 49);

    LeftRightPlanarityTest test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TEST_F(LeftRightPlanarityTestGTest, RandomGeneratedPlanarGraph100Nodes) {
    Graph graph(100);

    graph.addEdge(0, 92);
    graph.addEdge(1, 40);
    graph.addEdge(1, 49);
    graph.addEdge(1, 91);
    graph.addEdge(2, 87);
    graph.addEdge(3, 79);
    graph.addEdge(3, 84);
    graph.addEdge(3, 90);
    graph.addEdge(4, 50);
    graph.addEdge(5, 11);
    graph.addEdge(6, 40);
    graph.addEdge(6, 88);
    graph.addEdge(7, 48);
    graph.addEdge(7, 71);
    graph.addEdge(7, 92);
    graph.addEdge(8, 71);
    graph.addEdge(9, 14);
    graph.addEdge(9, 29);
    graph.addEdge(9, 91);
    graph.addEdge(10, 65);
    graph.addEdge(10, 93);
    graph.addEdge(12, 77);
    graph.addEdge(12, 94);
    graph.addEdge(13, 16);
    graph.addEdge(13, 28);
    graph.addEdge(14, 36);
    graph.addEdge(15, 42);
    graph.addEdge(16, 51);
    graph.addEdge(17, 75);
    graph.addEdge(19, 27);
    graph.addEdge(19, 33);
    graph.addEdge(19, 52);
    graph.addEdge(19, 60);
    graph.addEdge(19, 64);
    graph.addEdge(19, 91);
    graph.addEdge(20, 56);
    graph.addEdge(20, 84);
    graph.addEdge(21, 78);
    graph.addEdge(21, 89);
    graph.addEdge(22, 50);
    graph.addEdge(22, 51);
    graph.addEdge(22, 53);
    graph.addEdge(22, 75);
    graph.addEdge(23, 38);
    graph.addEdge(24, 65);
    graph.addEdge(25, 76);
    graph.addEdge(25, 90);
    graph.addEdge(26, 50);
    graph.addEdge(27, 35);
    graph.addEdge(27, 56);
    graph.addEdge(27, 84);
    graph.addEdge(28, 75);
    graph.addEdge(30, 34);
    graph.addEdge(30, 38);
    graph.addEdge(30, 51);
    graph.addEdge(31, 60);
    graph.addEdge(32, 35);
    graph.addEdge(32, 71);
    graph.addEdge(32, 76);
    graph.addEdge(33, 53);
    graph.addEdge(35, 95);
    graph.addEdge(35, 96);
    graph.addEdge(36, 53);
    graph.addEdge(36, 88);
    graph.addEdge(37, 68);
    graph.addEdge(39, 49);
    graph.addEdge(41, 55);
    graph.addEdge(41, 74);
    graph.addEdge(41, 75);
    graph.addEdge(42, 69);
    graph.addEdge(43, 46);
    graph.addEdge(43, 77);
    graph.addEdge(44, 77);
    graph.addEdge(45, 56);
    graph.addEdge(46, 67);
    graph.addEdge(47, 69);
    graph.addEdge(47, 79);
    graph.addEdge(47, 86);
    graph.addEdge(51, 63);
    graph.addEdge(51, 68);
    graph.addEdge(51, 72);
    graph.addEdge(51, 74);
    graph.addEdge(52, 75);
    graph.addEdge(52, 86);
    graph.addEdge(54, 72);
    graph.addEdge(58, 62);
    graph.addEdge(58, 63);
    graph.addEdge(58, 65);
    graph.addEdge(58, 66);
    graph.addEdge(58, 78);
    graph.addEdge(60, 87);
    graph.addEdge(61, 82);
    graph.addEdge(61, 95);
    graph.addEdge(65, 83);
    graph.addEdge(67, 82);
    graph.addEdge(68, 87);
    graph.addEdge(70, 80);
    graph.addEdge(72, 73);
    graph.addEdge(75, 80);
    graph.addEdge(77, 78);
    graph.addEdge(80, 89);
    graph.addEdge(84, 86);
    graph.addEdge(89, 95);
    LeftRightPlanarityTest test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TEST_F(LeftRightPlanarityTestGTest, RandomGeneratedPlanarConnectedGraph51Nodes) {
    Graph graph(51);
    graph.addEdge(0, 14);
    graph.addEdge(0, 17);
    graph.addEdge(0, 24);
    graph.addEdge(0, 43);
    graph.addEdge(0, 49);
    graph.addEdge(1, 15);
    graph.addEdge(1, 25);
    graph.addEdge(2, 38);
    graph.addEdge(2, 46);
    graph.addEdge(3, 8);
    graph.addEdge(3, 12);
    graph.addEdge(3, 22);
    graph.addEdge(4, 24);
    graph.addEdge(5, 44);
    graph.addEdge(5, 45);
    graph.addEdge(6, 36);
    graph.addEdge(6, 45);
    graph.addEdge(7, 9);
    graph.addEdge(7, 18);
    graph.addEdge(7, 20);
    graph.addEdge(7, 28);
    graph.addEdge(8, 44);
    graph.addEdge(10, 32);
    graph.addEdge(11, 40);
    graph.addEdge(11, 50);
    graph.addEdge(12, 47);
    graph.addEdge(13, 19);
    graph.addEdge(14, 41);
    graph.addEdge(14, 45);
    graph.addEdge(14, 49);
    graph.addEdge(15, 47);
    graph.addEdge(16, 19);
    graph.addEdge(17, 26);
    graph.addEdge(18, 30);
    graph.addEdge(19, 29);
    graph.addEdge(20, 47);
    graph.addEdge(21, 23);
    graph.addEdge(22, 29);
    graph.addEdge(22, 50);
    graph.addEdge(23, 28);
    graph.addEdge(24, 38);
    graph.addEdge(25, 28);
    graph.addEdge(25, 31);
    graph.addEdge(25, 32);
    graph.addEdge(26, 48);
    graph.addEdge(27, 29);
    graph.addEdge(27, 36);
    graph.addEdge(27, 50);
    graph.addEdge(28, 30);
    graph.addEdge(30, 32);
    graph.addEdge(31, 47);
    graph.addEdge(33, 36);
    graph.addEdge(33, 37);
    graph.addEdge(33, 48);
    graph.addEdge(34, 50);
    graph.addEdge(35, 46);
    graph.addEdge(37, 39);
    graph.addEdge(37, 42);
    graph.addEdge(38, 46);
    graph.addEdge(43, 50);
    graph.addEdge(44, 49);
    LeftRightPlanarityTest test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

TEST_F(LeftRightPlanarityTestGTest, TwoSubGraphsOnePlanarOneNonPlanar) {
    Graph graph(10);
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);
    graph.addEdge(0, 3);
    graph.addEdge(0, 4);
    graph.addEdge(0, 5);
    graph.addEdge(1, 3);
    graph.addEdge(1, 4);
    graph.addEdge(1, 5);
    graph.addEdge(2, 3);
    graph.addEdge(2, 4);
    graph.addEdge(2, 5);
    graph.addEdge(6, 7);
    graph.addEdge(7, 8);
    graph.addEdge(8, 9);
    graph.addEdge(9, 6);
    LeftRightPlanarityTest test(graph);
    test.run();
    EXPECT_FALSE(test.isPlanar());
}

TEST_F(LeftRightPlanarityTestGTest, TwoRandomSubGraphsBothPlanar) {
    Graph graph(20);
    graph.addEdge(0, 9);
    graph.addEdge(0, 12);
    graph.addEdge(0, 16);
    graph.addEdge(0, 17);
    graph.addEdge(1, 11);
    graph.addEdge(2, 6);
    graph.addEdge(2, 7);
    graph.addEdge(2, 19);
    graph.addEdge(3, 7);
    graph.addEdge(3, 15);
    graph.addEdge(4, 9);
    graph.addEdge(5, 18);
    graph.addEdge(6, 8);
    graph.addEdge(6, 10);
    graph.addEdge(6, 14);
    graph.addEdge(6, 19);
    graph.addEdge(10, 13);
    graph.addEdge(11, 18);
    graph.addEdge(16, 17);
    graph.addEdge(16, 18);
    graph.addEdge(17, 18);
    LeftRightPlanarityTest test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

} // namespace NetworKit
