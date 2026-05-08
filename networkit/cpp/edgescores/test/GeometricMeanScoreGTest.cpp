/*
 * GeometricMeanScoreGTest.cpp
 *
 * Created on: 08.05.2026
 * Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 */

#include <gtest/gtest.h>

#include <cmath>
#include <networkit/edgescores/GeometricMeanScore.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class GeometricMeanScoreGTest : public ::testing::Test {};

TEST_F(GeometricMeanScoreGTest, testRunThrowsIfEdgesAreNotIndexed) {
    Graph G(2, false, false);
    G.addEdge(0, 1);

    const std::vector<double> attribute{1.0};

    GeometricMeanScore score(G, attribute);

    EXPECT_THROW(score.run(), std::runtime_error);
}

TEST_F(GeometricMeanScoreGTest, testComputesScoresOnPath) {
    Graph G(4, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.indexEdges();

    std::vector<double> attribute(G.upperEdgeIdBound());
    attribute[G.edgeId(0, 1)] = 2.0;
    attribute[G.edgeId(1, 2)] = 4.0;
    attribute[G.edgeId(1, 3)] = 6.0;

    GeometricMeanScore score(G, attribute);
    score.run();

    const auto &scores = score.scores();

    EXPECT_DOUBLE_EQ(scores[G.edgeId(0, 1)], 2.0 / std::sqrt(2.0 * 12.0));
    EXPECT_DOUBLE_EQ(scores[G.edgeId(1, 2)], 4.0 / std::sqrt(12.0 * 4.0));
    EXPECT_DOUBLE_EQ(scores[G.edgeId(1, 3)], 6.0 / std::sqrt(12.0 * 6.0));
}

TEST_F(GeometricMeanScoreGTest, testKeepsZeroAttributeScoreAtZero) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.indexEdges();

    std::vector<double> attribute(G.upperEdgeIdBound());
    attribute[G.edgeId(0, 1)] = 0.0;
    attribute[G.edgeId(1, 2)] = 4.0;

    GeometricMeanScore score(G, attribute);
    score.run();

    const auto &scores = score.scores();

    EXPECT_DOUBLE_EQ(scores[G.edgeId(0, 1)], 0.0);
    EXPECT_DOUBLE_EQ(scores[G.edgeId(1, 2)], 1.0);
}

} // namespace NetworKit
