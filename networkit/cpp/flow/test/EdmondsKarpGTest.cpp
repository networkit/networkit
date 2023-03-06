/*
 * EdmondsKarpGTest.cpp
 *
 *  Created on: 13.06.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <networkit/flow/EdmondsKarp.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class EdmondsKarpGTest : public testing::Test {};

TEST_F(EdmondsKarpGTest, testEdmondsKarpP1) {
    Graph G(7, false);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(0, 3);
    G.addEdge(1, 2);
    G.addEdge(1, 4);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 6);
    G.addEdge(5, 6);

    G.indexEdges();

    EdmondsKarp edKa(G, 0, 6);
    edKa.run();
    EXPECT_DOUBLE_EQ(2, edKa.getMaxFlow()) << "max flow is not correct";

    EXPECT_DOUBLE_EQ(1, edKa.getFlow(4, 6));
    EXPECT_DOUBLE_EQ(1, edKa.getFlow(5, 6));

    std::vector<node> sourceSet(edKa.getSourceSet());
    EXPECT_THAT(sourceSet, testing::UnorderedElementsAre(0, 1, 2, 3, 4));
    EXPECT_THAT(sourceSet, testing::Not(testing::Contains(testing::AnyOf(5, 6))));
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpP2) {
    Graph G(6, true);
    G.addEdge(0, 1, 5);
    G.addEdge(0, 2, 15);
    G.addEdge(1, 3, 5);
    G.addEdge(1, 4, 5);
    G.addEdge(2, 3, 5);
    G.addEdge(2, 4, 5);
    G.addEdge(3, 5, 15);
    G.addEdge(4, 5, 5);

    G.indexEdges();

    EdmondsKarp edKa(G, 0, 5);
    edKa.run();

    EXPECT_DOUBLE_EQ(15, edKa.getMaxFlow()) << "max flow is not correct";
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpP3) {
    Graph G(5, false, false);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(1, 3);
    G.addEdge(2, 3);
    G.addEdge(3, 4);
    G.indexEdges();

    EdmondsKarp edKa(G, 0, 4);
    edKa.run();
    EXPECT_DOUBLE_EQ(1, edKa.getMaxFlow()) << "max flow is not correct";

    std::vector<node> sourceSet(edKa.getSourceSet());
    EXPECT_THAT(sourceSet, testing::UnorderedElementsAre(0, 1, 2, 3));
    EXPECT_THAT(sourceSet, testing::Not(::testing::Contains(4)));
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpDirected1) {
    Graph G(4, true, true);
    G.addEdge(0, 1, 10);
    G.addEdge(0, 2, 5);
    G.addEdge(1, 2, 15);
    G.addEdge(1, 3, 5);
    G.addEdge(2, 3, 10);
    G.indexEdges();

    EdmondsKarp edKa(G, 0, 3);
    edKa.run();
    EXPECT_DOUBLE_EQ(15, edKa.getMaxFlow()) << "max flow is not correct";

    std::vector<node> sourceSet(edKa.getSourceSet());
    EXPECT_THAT(sourceSet, testing::Contains(0));
    // either 1 and 2 are in the set or none
    EXPECT_EQ(std::count(sourceSet.begin(), sourceSet.end(), 1),
              std::count(sourceSet.begin(), sourceSet.end(), 2));
    EXPECT_THAT(sourceSet, testing::Not(::testing::Contains(3)));
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpDirected2) {
    Graph G(5, false, true);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(1, 3);
    G.addEdge(2, 3);
    G.addEdge(3, 4);
    G.indexEdges();

    EdmondsKarp edKa(G, 0, 4);
    edKa.run();
    EXPECT_DOUBLE_EQ(1, edKa.getMaxFlow()) << "max flow is not correct";

    std::vector<node> sourceSet(edKa.getSourceSet());
    EXPECT_THAT(sourceSet, testing::UnorderedElementsAre(0, 1, 2, 3));
    EXPECT_THAT(sourceSet, testing::Not(testing::Contains(4)));
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpDirected3) {
    Graph G(9, true, true);
    G.addEdge(0, 1, 2);
    G.addEdge(0, 2, 2);
    G.addEdge(0, 3, 2);
    G.addEdge(1, 4, 2);
    G.addEdge(2, 4, 2);
    G.addEdge(3, 4, 2);
    G.addEdge(4, 5, 4);
    G.addEdge(5, 6, 2);
    G.addEdge(5, 7, 2);
    G.addEdge(6, 8, 2);
    G.addEdge(7, 8, 2);
    G.indexEdges();

    EdmondsKarp edKa(G, 0, 8);
    edKa.run();
    EXPECT_DOUBLE_EQ(4, edKa.getMaxFlow()) << "max flow is not correct";

    std::vector<node> sourceSet(edKa.getSourceSet());
    EXPECT_THAT(sourceSet, testing::UnorderedElementsAre(0, 1, 2, 3, 4));
    EXPECT_THAT(sourceSet, testing::Not(testing::Contains(8)));
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpDirected4) {
    Graph G(9, true, true);
    G.addEdge(8, 7, 2);
    G.addEdge(8, 6, 2);
    G.addEdge(8, 5, 2);
    G.addEdge(7, 4, 2);
    G.addEdge(6, 4, 2);
    G.addEdge(5, 4, 2);
    G.addEdge(4, 3, 4);
    G.addEdge(3, 2, 2);
    G.addEdge(3, 1, 2);
    G.addEdge(2, 0, 2);
    G.addEdge(1, 0, 2);
    G.indexEdges();

    EdmondsKarp edKa(G, 8, 0);
    edKa.run();
    EXPECT_DOUBLE_EQ(4, edKa.getMaxFlow()) << "max flow is not correct";

    std::vector<node> sourceSet(edKa.getSourceSet());
    EXPECT_THAT(sourceSet, testing::Contains(8));
    EXPECT_THAT(sourceSet, testing::Not(testing::Contains(0)));
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpDirected5) {
    Graph G(4, true, true);
    G.addEdge(3, 2, 3.4);
    G.addEdge(2, 1, 2.4);
    G.addEdge(1, 0, 4.4);
    G.indexEdges();

    EdmondsKarp edKa(G, 3, 0);
    edKa.run();
    EXPECT_DOUBLE_EQ(2.4, edKa.getMaxFlow()) << "max flow is not correct";

    std::vector<node> sourceSet(edKa.getSourceSet());
    EXPECT_THAT(sourceSet, testing::Contains(3));
    EXPECT_THAT(sourceSet, testing::Not(testing::Contains(0)));
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpDirected6) {
    Graph G(4, true, true);
    G.addEdge(3, 2, 3.4);
    G.addEdge(1, 2, 2.4);
    G.addEdge(1, 0, 4.4);
    G.indexEdges();

    EdmondsKarp edKa(G, 3, 0);
    edKa.run();
    EXPECT_DOUBLE_EQ(0, edKa.getMaxFlow()) << "max flow is not correct";

    std::vector<node> sourceSet(edKa.getSourceSet());
    EXPECT_THAT(sourceSet, testing::Contains(3));
    EXPECT_THAT(sourceSet, testing::Not(testing::Contains(0)));
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpDirected7) {
    // Example from wikipedia
    // (https://en.wikipedia.org/wiki/Edmonds%E2%80%93Karp_algorithm#Example)
    Graph G(7, true, true);
    G.addEdge(0, 1, 3);
    G.addEdge(0, 3, 3);
    G.addEdge(1, 2, 4);
    G.addEdge(2, 0, 3);
    G.addEdge(2, 3, 1);
    G.addEdge(2, 4, 2);
    G.addEdge(3, 4, 2);
    G.addEdge(3, 5, 6);
    G.addEdge(4, 1, 1);
    G.addEdge(4, 6, 1);
    G.addEdge(5, 6, 9);
    G.indexEdges();

    EdmondsKarp edKa(G, 0, 6);
    edKa.run();
    EXPECT_DOUBLE_EQ(5, edKa.getMaxFlow()) << "max flow is not correct";

    std::vector<node> sourceSet(edKa.getSourceSet());
    EXPECT_THAT(sourceSet, testing::Contains(0));
    EXPECT_THAT(sourceSet, testing::Not(testing::Contains(6)));
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpTwoPaths) {
    Graph G(11);

    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(2, 3);
    G.addEdge(3, 4);
    G.addEdge(1, 4);
    G.addEdge(1, 5);
    G.addEdge(5, 6);
    G.addEdge(6, 7);
    G.addEdge(7, 8);
    G.addEdge(8, 9);
    G.addEdge(4, 10);
    G.addEdge(9, 10);

    G.indexEdges();

    EdmondsKarp edKa(G, 0, 10);
    edKa.run();

    EXPECT_DOUBLE_EQ(2, edKa.getMaxFlow());
    EXPECT_DOUBLE_EQ(0, edKa.getFlow(1, 4));
}

TEST_F(EdmondsKarpGTest, testEdmondsKarpUnconnected) {
    Graph G(6, true);
    G.addEdge(0, 1, 5);
    G.addEdge(0, 2, 15);
    G.addEdge(1, 2, 5);
    G.addEdge(3, 4, 5);
    G.addEdge(3, 5, 15);
    G.addEdge(4, 5, 5);

    G.indexEdges();

    EdmondsKarp edKa(G, 0, 5);
    edKa.run();
    EXPECT_DOUBLE_EQ(0, edKa.getMaxFlow()) << "max flow is not correct";
}

} /* namespace NetworKit */
