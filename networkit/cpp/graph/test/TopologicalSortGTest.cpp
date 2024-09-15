/*
 * TopologicalSort.hpp
 *
 *  Created on: 22.11.2021
 *      Author: Fabian Brandt-Tumescheit
 */

#include <networkit/graph/TopologicalSort.hpp>

#include <stdexcept>
#include <gtest/gtest.h>

namespace NetworKit {

class TopologicalSortGTest : public testing::Test {
protected:
    Graph inputGraph(bool directed) const noexcept;
};

Graph TopologicalSortGTest::inputGraph(bool directed) const noexcept {
    auto G = Graph(5, false, directed);

    /**
     * /--> 1 --> 3
     * |    ^
     * 0    |
     * |    |
     * \--> 2 <-- 4
     */

    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(2, 1);
    G.addEdge(1, 3);
    G.addEdge(4, 2);

    return G;
}

TEST_F(TopologicalSortGTest, testTopologicalSort) {
    auto G = inputGraph(true);

    TopologicalSort topSort = TopologicalSort(G);
    topSort.run();
    std::vector<node> res = topSort.getResult();

    // Test for valid topology
    auto indexNode0 = std::distance(res.begin(), std::find(res.begin(), res.end(), 0));
    auto indexNode1 = std::distance(res.begin(), std::find(res.begin(), res.end(), 1));
    auto indexNode2 = std::distance(res.begin(), std::find(res.begin(), res.end(), 2));
    auto indexNode3 = std::distance(res.begin(), std::find(res.begin(), res.end(), 3));
    auto indexNode4 = std::distance(res.begin(), std::find(res.begin(), res.end(), 4));
    EXPECT_EQ(res.size(), 5);
    // node 2 is depending on node 0
    EXPECT_TRUE(indexNode2 > indexNode0);
    // node 2 is depending on node 4
    EXPECT_TRUE(indexNode2 > indexNode4);
    // node 1 is depending on node 2
    EXPECT_TRUE(indexNode1 > indexNode2);
    // node 3 is depending on node 1
    EXPECT_TRUE(indexNode3 > indexNode1);
}

TEST_F(TopologicalSortGTest, testRepeatedRuns) {
    auto G = inputGraph(true);
    TopologicalSort topSort = TopologicalSort(G);
    topSort.run();
    std::vector<node> res = topSort.getResult();
    topSort.run();
    std::vector<node> res2 = topSort.getResult();
    EXPECT_TRUE(res == res2);
}

TEST_F(TopologicalSortGTest, testRejectGraphWithCycles) {
    auto G = inputGraph(true);
    G.addEdge(3, 4);
    TopologicalSort topSort = TopologicalSort(G);
    EXPECT_THROW(topSort.run(), std::runtime_error);
}

TEST_F(TopologicalSortGTest, testRejectUndirectedGraph) {
    EXPECT_THROW(TopologicalSort(inputGraph(false)), std::runtime_error);
}

TEST_F(TopologicalSortGTest, testNonContinuousNodeIds) {
    auto G = inputGraph(true);
    G.removeNode(3);

    TopologicalSort topSort = TopologicalSort(G);
    topSort.run();
    std::vector<node> res = topSort.getResult();

    // Test for valid topology
    auto indexNode0 = std::distance(res.begin(), std::find(res.begin(), res.end(), 0));
    auto indexNode1 = std::distance(res.begin(), std::find(res.begin(), res.end(), 1));
    auto indexNode2 = std::distance(res.begin(), std::find(res.begin(), res.end(), 2));
    auto indexNode4 = std::distance(res.begin(), std::find(res.begin(), res.end(), 4));
    EXPECT_EQ(res.size(), 4);
    // node 2 is depending on node 0
    EXPECT_TRUE(indexNode2 > indexNode0);
    // node 2 is depending on node 4
    EXPECT_TRUE(indexNode2 > indexNode4);
    // node 1 is depending on node 2
    EXPECT_TRUE(indexNode1 > indexNode2);
}

TEST_F(TopologicalSortGTest, testCustomNodeIdMapping) {
    auto G = inputGraph(true);
    G.removeNode(3);
    std::unordered_map<node, node> mapping;
    mapping[0] = 0;
    mapping[1] = 1;
    mapping[2] = 2;
    mapping[4] = 3;
    TopologicalSort topSort = TopologicalSort(G, mapping);
    topSort.run();
    std::vector<node> res = topSort.getResult();

    // Test for valid topology
    auto indexNode0 = std::distance(res.begin(), std::find(res.begin(), res.end(), 0));
    auto indexNode1 = std::distance(res.begin(), std::find(res.begin(), res.end(), 1));
    auto indexNode2 = std::distance(res.begin(), std::find(res.begin(), res.end(), 2));
    auto indexNode4 = std::distance(res.begin(), std::find(res.begin(), res.end(), 4));
    EXPECT_EQ(res.size(), 4);
    // node 2 is depending on node 0
    EXPECT_TRUE(indexNode2 > indexNode0);
    // node 2 is depending on node 4
    EXPECT_TRUE(indexNode2 > indexNode4);
    // node 1 is depending on node 2
    EXPECT_TRUE(indexNode1 > indexNode2);
}

TEST_F(TopologicalSortGTest, testWrongSizeOfNodeIdMapping) {
    auto G = inputGraph(true);
    G.removeNode(3);
    std::unordered_map<node, node> mapping;
    mapping[0] = 0;
    mapping[1] = 1;
    mapping[2] = 2;
    mapping[4] = 3;
    mapping[5] = 4;
    EXPECT_THROW(TopologicalSort(G, mapping), std::runtime_error);
}

TEST_F(TopologicalSortGTest, testNonContinuousNodeIdMapping) {
    auto G = inputGraph(true);
    G.removeNode(3);
    std::unordered_map<node, node> mapping;
    mapping[0] = 0;
    mapping[1] = 4;
    mapping[2] = 2;
    mapping[4] = 3;
    EXPECT_THROW(TopologicalSort(G, mapping, true), std::runtime_error);
}

TEST_F(TopologicalSortGTest, testNonInjectiveNodeIdMapping) {
    auto G = inputGraph(true);
    G.removeNode(3);
    std::unordered_map<node, node> mapping;
    mapping[0] = 0;
    mapping[1] = 1;
    mapping[2] = 1;
    mapping[4] = 3;
    EXPECT_THROW(TopologicalSort(G, mapping, true), std::runtime_error);
}
} // namespace NetworKit
