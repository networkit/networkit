//
// Created by Alexander K. Ziebs on 25.07.2024.
//

#include <gtest/gtest.h>
#include <networkit/embedding/KHop.hpp>

namespace NetworKit {

class KHopTest : public testing::Test {};

TEST_F(KHopTest, testConstructStrictKHopGraph) {

    // pre-Hop Graph:
    /*
     *        2 --- 1
     *        |   / |
     *        |  /  |
     *        | /   |
     *  4 --- 3     0
     */

    Graph preHop(5, false, false);
    preHop.addEdge(0, 1);
    preHop.addEdge(1, 2);
    preHop.addEdge(1, 3);
    preHop.addEdge(2, 3);
    preHop.addEdge(3, 4);

    KHop hop(preHop, 2, 6.25, 80, 10, 128, KHop::khopMode::STRICT);
    Graph G_2 = hop.GetG_k();

    // Expected 2 Hop Graph G_2:
    /*
     *  4 --- 2 --- 0
     *  |           |
     *  |           |
     *  1           3
     */

    // Graph expected(5, false, false);
    // G.addEdge(0, 3);
    // G.addEdge(0, 2);
    // G.addEdge(2, 4);
    // G.addEdge(4, 1);

    // check for all edges, that should exist
    EXPECT_TRUE(G_2.hasEdge(0, 3));
    EXPECT_TRUE(G_2.hasEdge(0, 2));
    EXPECT_TRUE(G_2.hasEdge(2, 4));
    EXPECT_TRUE(G_2.hasEdge(4, 1));

    // check that no additional edges exist
    EXPECT_TRUE(G_2.numberOfEdges() == 4);
}

TEST_F(KHopTest, testConstructDefaultKHopGraph) {

    // pre-Hop Graph:
    /*
     *  2 --- 1
     *      / |
     *     /  |
     *    /   |
     *  3     0
     */

    Graph preHop(4, false, false);
    preHop.addEdge(0, 1);
    preHop.addEdge(1, 2);
    preHop.addEdge(1, 3);

    KHop hop(preHop, 2, 6.25, 80, 10, 128, KHop::khopMode::DEFAULT);
    Graph G_2 = hop.GetG_k();

    // Expected 2 Hop Graph G_2 should be a fully connected Graph with 4 nodes (including self
    // loops)

    // check for all edges, that should exist
    EXPECT_TRUE(G_2.hasEdge(0, 0));
    EXPECT_TRUE(G_2.hasEdge(0, 1));
    EXPECT_TRUE(G_2.hasEdge(0, 2));
    EXPECT_TRUE(G_2.hasEdge(0, 3));
    EXPECT_TRUE(G_2.hasEdge(1, 0));
    EXPECT_TRUE(G_2.hasEdge(1, 1));
    EXPECT_TRUE(G_2.hasEdge(1, 2));
    EXPECT_TRUE(G_2.hasEdge(1, 3));
    EXPECT_TRUE(G_2.hasEdge(2, 0));
    EXPECT_TRUE(G_2.hasEdge(2, 1));
    EXPECT_TRUE(G_2.hasEdge(2, 2));
    EXPECT_TRUE(G_2.hasEdge(2, 3));
    EXPECT_TRUE(G_2.hasEdge(3, 0));
    EXPECT_TRUE(G_2.hasEdge(3, 1));
    EXPECT_TRUE(G_2.hasEdge(3, 2));
    EXPECT_TRUE(G_2.hasEdge(3, 3));
}

} // namespace NetworKit
