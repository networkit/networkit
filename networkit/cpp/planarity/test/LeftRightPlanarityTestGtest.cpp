//
// Created by andreas on 05.01.25.
//
#include <gtest/gtest.h>
#include <networkit/planarity/LeftRightPlanarityTest.hpp>

namespace NetworKit {
class LeftRightPlanarityTestGTest : public testing::Test {
public:
    LeftRightPlanarityTestGTest() = default;
    static constexpr int max_number_of_nodes{50};
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

};

TEST_F(LeftRightPlanarityTestGTest, PathGraph) {

    for (count numberOfNodes{2}; numberOfNodes <= max_number_of_nodes; ++numberOfNodes)
    {

        auto graph = pathGraph(numberOfNodes);

        LeftRightPlanarityTest test(graph);
        test.run();
        EXPECT_TRUE(test.isPlanar());
    }

}
}