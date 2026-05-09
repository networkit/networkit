/*
 * PrefixJaccardScoreGTest.cpp
 *
 * Created on: 08.05.2026
 * Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/edgescores/PrefixJaccardScore.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

using testing::DoubleEq;
using testing::ElementsAre;

class PrefixJaccardScoreGTest : public testing::Test {};

TEST_F(PrefixJaccardScoreGTest, testRunThrowsIfEdgesAreNotIndexed) {
    Graph G(2, false, false);
    G.addEdge(0, 1);

    const std::vector<double> attribute{1.0};

    PrefixJaccardScore<double> score(G, attribute);

    EXPECT_THROW(score.run(), std::runtime_error);
}

TEST_F(PrefixJaccardScoreGTest, testTriangleHasFullPrefixOverlap) {
    Graph G(4, false, false);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.indexEdges();

    std::vector<double> attribute(G.upperEdgeIdBound(), 1.0);

    PrefixJaccardScore<double> score(G, attribute);
    score.run();

    EXPECT_THAT(score.scores(), ElementsAre(DoubleEq(1.0), DoubleEq(1.0), DoubleEq(1.0)));
}

TEST_F(PrefixJaccardScoreGTest, testPathHasNoPrefixOverlap) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.indexEdges();

    std::vector<double> attribute(G.upperEdgeIdBound(), 1.0);

    PrefixJaccardScore<double> score(G, attribute);
    score.run();

    EXPECT_THAT(score.scores(), ElementsAre(DoubleEq(0.0), DoubleEq(0.0)));
}

TEST_F(PrefixJaccardScoreGTest, testBestPrefixCanBeBetterThanFullNeighborhoodOverlap) {
    Graph G(5, false, false);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(0, 3);
    G.addEdge(1, 4);
    G.indexEdges();

    std::vector<double> attribute(G.upperEdgeIdBound());
    attribute[G.edgeId(0, 1)] = 5.0;
    attribute[G.edgeId(0, 2)] = 10.0;
    attribute[G.edgeId(1, 2)] = 10.0;
    attribute[G.edgeId(0, 3)] = 1.0;
    attribute[G.edgeId(1, 4)] = 1.0;

    PrefixJaccardScore<double> score(G, attribute);
    score.run();

    EXPECT_THAT(score.scores(), ElementsAre(DoubleEq(1.0), DoubleEq(1.0), DoubleEq(1.0),
                                            DoubleEq(0.0), DoubleEq(0.0)));
}

TEST_F(PrefixJaccardScoreGTest, testComputesFractionalPrefixOverlap) {
    Graph G(4, false, false);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.indexEdges();

    std::vector<double> attribute(G.upperEdgeIdBound(), 1.0);

    PrefixJaccardScore<double> score(G, attribute);
    score.run();

    EXPECT_THAT(score.scores(),
                ElementsAre(DoubleEq(0.5), DoubleEq(1.0), DoubleEq(0.5), DoubleEq(0.0)));
}
} // namespace NetworKit
