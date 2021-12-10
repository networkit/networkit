/*
 * TopologicalSort.hpp
 *
 *  Created on: 22.11.2021
 *      Author: Fabian Brandt-Tumescheit
 */

#include <networkit/graph/TopologicalSort.hpp>

#include <iostream>
#include <gtest/gtest.h>

namespace NetworKit {

class TopologicalSortGTest : public testing::TestWithParam<bool> {
protected:
    bool directed() const noexcept;
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, TopologicalSortGTest, testing::Values(false, true));

bool TopologicalSortGTest::directed() const noexcept {
    return GetParam();
}

TEST_P(TopologicalSortGTest, testTopologicalSort) {
    auto G = Graph(5, false, directed());

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

    if (directed()) {
        TopologicalSort topSort = TopologicalSort(G);
        topSort.run();
        std::vector<node> res = topSort.getResult();

        // Test for valid topology
        auto indexNode0 = std::distance(res.begin(), std::find(res.begin(), res.end(), 0));
        auto indexNode1 = std::distance(res.begin(), std::find(res.begin(), res.end(), 1));
        auto indexNode2 = std::distance(res.begin(), std::find(res.begin(), res.end(), 2));
        auto indexNode3 = std::distance(res.begin(), std::find(res.begin(), res.end(), 3));
        auto indexNode4 = std::distance(res.begin(), std::find(res.begin(), res.end(), 4));
        // node 2 is depending on node 0
        EXPECT_TRUE(indexNode2 > indexNode0);
        // node 2 is depending on node 4
        EXPECT_TRUE(indexNode2 > indexNode4);
        // node 1 is depending on node 2
        EXPECT_TRUE(indexNode1 > indexNode2);
        // node 3 is depending on node 1
        EXPECT_TRUE(indexNode3 > indexNode1);

        // Test for repeated runs
        topSort.run();
        std::vector<node> res2 = topSort.getResult();
        EXPECT_TRUE(res == res2);

        // Test cycle detection
        G.addEdge(3, 4);
        topSort = TopologicalSort(G);
        EXPECT_ANY_THROW(topSort.run());
    } else {
        EXPECT_ANY_THROW(new TopologicalSort(G));
    }
}
} // namespace NetworKit