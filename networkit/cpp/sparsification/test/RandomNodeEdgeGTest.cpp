/*
 * RandomNodeEdgeGTest.cpp
 *
 *  Created on: 05.09.2023
 *      Author: Angus Petholz
 */

#include <gtest/gtest.h>

#include <networkit/sparsification/RandomEdgeScore.hpp>
#include <networkit/sparsification/RandomNodeEdgeScore.hpp>

namespace NetworKit {

class RandomNodeEdgeGTest : public testing::Test {};

Graph initGraphUnindexed() {
    Graph G(5);      // EdgeId
    G.addEdge(0, 1); // 0     G:  0 - 1
    G.addEdge(0, 2); // 1         | / |
    G.addEdge(1, 2); // 2         2 - 3
    G.addEdge(1, 3); // 3         |
    G.addEdge(2, 3); // 4         4
    G.addEdge(2, 4); // 5
    return G;
}

TEST_F(RandomNodeEdgeGTest, testRandomEdgeScores) {
    Graph G = initGraphUnindexed();
    G.indexEdges();

    RandomEdgeScore rndes(G);
    rndes.run();
    auto res = rndes.scores();
    EXPECT_EQ(res.size(), G.numberOfEdges());
    for (auto it : res) {
        EXPECT_NEAR(it, 0.5, 0.5); // all entries in range[0:1]
    }
}

TEST_F(RandomNodeEdgeGTest, testRandomEdgeScore) {
    Graph G = initGraphUnindexed();
    G.indexEdges();

    RandomEdgeScore rndes(G);
    rndes.run();
    auto resNode = rndes.score(0, 1);
    auto resEid = rndes.score(0);

    EXPECT_EQ(resNode, resEid);
    EXPECT_NEAR(resNode, 0.5, 0.5); // entry in range[0:1]
}

TEST_F(RandomNodeEdgeGTest, testRandomEdgeScoreUnindexedEdges) {
    Graph G = initGraphUnindexed();

    RandomEdgeScore rndes(G);
    EXPECT_THROW(rndes.run(), std::runtime_error);
}

TEST_F(RandomNodeEdgeGTest, testRandomNodeEdgeScores) {
    Graph G = initGraphUnindexed();
    G.indexEdges();

    RandomNodeEdgeScore rndnes(G);
    rndnes.run();
    auto res = rndnes.scores();
    EXPECT_EQ(res.size(), G.numberOfEdges());
    for (auto it : res) {
        EXPECT_NEAR(it, 0.5, 0.5); // all entries in range[0:1]
    }
}

TEST_F(RandomNodeEdgeGTest, testRandomNodeEdgeScore) {
    Graph G = initGraphUnindexed();
    G.indexEdges();

    RandomNodeEdgeScore rndnes(G);
    rndnes.run();
    auto resNode = rndnes.score(0, 1);
    auto resEid = rndnes.score(0);

    EXPECT_EQ(resNode, resEid);
    EXPECT_NEAR(resNode, 0.5, 0.5); // entry in range[0:1]
}

TEST_F(RandomNodeEdgeGTest, testRandomNodeEdgeScoreUnindexedEdges) {
    Graph G = initGraphUnindexed();

    RandomEdgeScore rndes(G);
    EXPECT_THROW(rndes.run(), std::runtime_error);
}

} // namespace NetworKit
/* namespace NetworKit */
