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
    static Graph pathGraph(const count numNodes)
    {
        Graph graph;
        graph.addNodes(numNodes);
        for (count i{}; i < numNodes-1; ++i)
        {
            graph.addEdge(i, i+1);
        }
        return graph;
    }

    static Graph cycleGraph(const count numNodes)
    {
        Graph graph;
        graph.addNodes(numNodes);
        for (count i{}; i < numNodes-1; ++i)
        {
            graph.addEdge(i, i+1);
        }
        if (numNodes > 2)
            graph.addEdge(numNodes-2, 0);
        return graph;
    }

    static Graph starGraph(const count numNodes)
    {
        Graph graph;
        graph.addNodes(numNodes);
        for (count i{}; i < numNodes-1; ++i)
        {
            graph.addEdge(i, i+1);
        }
        if (numNodes > 2)
            graph.addEdge(numNodes-2, 0);
        return graph;
    }

    static Graph binaryTreeGraph(const count numNodes)
    {
        Graph graph(numNodes);
        for (count i{}; i < numNodes; ++i)
        {
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

    static Graph wheelGraph(const count numNodes)
    {
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

    for (count numberOfNodes{2}; numberOfNodes <= max_number_of_nodes; ++numberOfNodes)
    {
        auto graph = pathGraph(numberOfNodes);
        LeftRightPlanarityTest test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityTestGTest, CycleGraph) {

    for (count numberOfNodes{2}; numberOfNodes <= max_number_of_nodes; ++numberOfNodes)
    {
        auto graph = cycleGraph(numberOfNodes);
        LeftRightPlanarityTest test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityTestGTest, StarGraph) {

    for (count numberOfNodes{2}; numberOfNodes <= max_number_of_nodes; ++numberOfNodes)
    {
        auto graph = starGraph(numberOfNodes);
        LeftRightPlanarityTest test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityTestGTest, TreeGraph) {

    for (count numberOfNodes{2}; numberOfNodes <= max_number_of_nodes; ++numberOfNodes)
    {
        auto graph = binaryTreeGraph(numberOfNodes);
        LeftRightPlanarityTest test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityTestGTest, WheelGraph) {

    for (count numberOfNodes{4}; numberOfNodes <= max_number_of_nodes; ++numberOfNodes)
    {
        auto graph = wheelGraph(numberOfNodes);
        LeftRightPlanarityTest test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }
}

TEST_F(LeftRightPlanarityTestGTest, GridGraph3x3) {

    Graph graph(9);

    // Horizontal edges
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(3, 4);
    graph.addEdge(4, 5);
    graph.addEdge(6, 7);
    graph.addEdge(7, 8);

    // Vertical edges
    graph.addEdge(0, 3);
    graph.addEdge(1, 4);
    graph.addEdge(2, 5);
    graph.addEdge(3, 6);
    graph.addEdge(4, 7);
    graph.addEdge(5, 8);
    LeftRightPlanarityTest test(graph);
    test.run();
    EXPECT_TRUE(test.isPlanar());
}

}