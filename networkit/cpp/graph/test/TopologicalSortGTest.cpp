/*
 * TopologicalSort.hpp
 *
 *  Created on: 22.11.2021
 *      Author: Fabian Brandt-Tumescheit
 */

#include <networkit/graph/TopologicalSort.hpp>

#include <stdexcept>
#include <unordered_map>
#include <gtest/gtest.h>

namespace NetworKit {

class TopologicalSortGTest : public testing::Test {
protected:
    Graph inputGraph(bool directed) const noexcept;
    std::unordered_map<node, node> makeMapping() const noexcept;
    void assertTopological(Graph &G, std::vector<node> &sort);
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

std::unordered_map<node, node> TopologicalSortGTest::makeMapping() const noexcept {
    std::unordered_map<node, node> mapping;
    mapping[0] = 0;
    mapping[1] = 1;
    mapping[2] = 2;
    mapping[4] = 3;

    return mapping;
}

void TopologicalSortGTest::assertTopological(Graph &G, std::vector<node> &sort) {
    std::unordered_map<node, node> indices;
    EXPECT_EQ(sort.size(), G.numberOfNodes());
    G.forNodes([&](node u) {
        indices[u] = std::distance(sort.begin(), std::find(sort.begin(), sort.end(), u));
    });
    G.forNodes([&](node u) {
        G.forNeighborsOf(u, [&](node v) { EXPECT_TRUE(indices[u] < indices[v]); });
    });
}

TEST_F(TopologicalSortGTest, testTopologicalSort) {
    auto G = inputGraph(true);

    TopologicalSort topSort = TopologicalSort(G);
    topSort.run();
    std::vector<node> res = topSort.getResult();

    assertTopological(G, res);
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

    assertTopological(G, res);
}

TEST_F(TopologicalSortGTest, testCustomNodeIdMapping) {
    auto G = inputGraph(true);
    G.removeNode(3);
    std::unordered_map<node, node> mapping = makeMapping();
    TopologicalSort topSort = TopologicalSort(G, mapping);
    topSort.run();
    std::vector<node> res = topSort.getResult();

    assertTopological(G, res);
}

TEST_F(TopologicalSortGTest, testWrongSizeOfNodeIdMapping) {
    auto G = inputGraph(true);
    G.removeNode(3);
    std::unordered_map<node, node> mapping = makeMapping();
    mapping[5] = 4;
    EXPECT_THROW(TopologicalSort(G, mapping), std::runtime_error);
}

TEST_F(TopologicalSortGTest, testNonContinuousNodeIdMapping) {
    auto G = inputGraph(true);
    G.removeNode(3);
    std::unordered_map<node, node> mapping = makeMapping();
    mapping[1] = 4;
    EXPECT_THROW(TopologicalSort(G, mapping, true), std::runtime_error);

    mapping = makeMapping();
    mapping.erase(1);
    // to get correct size
    mapping[5] = 5;
    EXPECT_THROW(TopologicalSort(G, mapping, true), std::runtime_error);
    TopologicalSort topSort = TopologicalSort(G, mapping);
    EXPECT_THROW(topSort.run(), std::runtime_error);
}

TEST_F(TopologicalSortGTest, testNonInjectiveNodeIdMapping) {
    auto G = inputGraph(true);
    G.removeNode(3);
    std::unordered_map<node, node> mapping = makeMapping();
    mapping[2] = 1;
    EXPECT_THROW(TopologicalSort(G, mapping, true), std::runtime_error);
}
} // namespace NetworKit
