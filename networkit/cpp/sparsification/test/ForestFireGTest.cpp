/*
 * ForestFireGTest.cpp
 *
 *  Created on: 05.09.2023
 *      Author: Angus Petholz
 */

#include <gtest/gtest.h>

#include <networkit/sparsification/ForestFireScore.hpp>

namespace NetworKit {

class ForestFireGTest : public testing::Test {};

TEST_F(ForestFireGTest, testForestFireScore) {
    Graph G(5);
    G.addEdge(0, 1); // 0     G:  0 - 1
    G.addEdge(0, 2); // 1         | / |
    G.addEdge(1, 2); // 2         2 - 3
    G.addEdge(1, 3); // 3         |
    G.addEdge(2, 3); // 4         4
    G.addEdge(2, 4); // 5
    G.indexEdges();  // EdgeId

    ForestFireScore ff(G, 0.5, 1);
    ff.run();
    auto res = ff.scores();
    EXPECT_EQ(res.size(), G.numberOfEdges());
}

} // namespace NetworKit
/* namespace NetworKit */
