/*
 * EdgeScoreGTest.cpp
 *
 *  Created on: 08.11.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 */

#include <gtest/gtest.h>

#include <networkit/edgescores/EdgeScore.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class EdgeScoreTest : public ::testing::Test {
protected:
    void compareTopology(const Graph &original, const Graph &result) {
        EXPECT_EQ(result.numberOfNodes(), original.numberOfNodes());
        EXPECT_EQ(result.numberOfEdges(), original.numberOfEdges());
        EXPECT_TRUE(result.isDirected() == original.isDirected());
    }
};

class TestEdgeScoreHelper final : public EdgeScore<double> {
public:
    explicit TestEdgeScoreHelper(const Graph &G) : EdgeScore<double>(G) {}

    void run() override {
        if (!G->hasEdgeIds()) {
            scoreData.clear();
        } else {
            scoreData.assign(G->upperEdgeIdBound(), 0.0);
            G->parallelForEdges(
                [&](node, node, edgeid eid) { scoreData[eid] = static_cast<double>(eid + 1); });
        }

        hasRun = true;
    }
};

TEST_F(EdgeScoreTest, testScoresAndScoreRequireRun) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.indexEdges();

    EdgeScore<double> baseScore(G);

    // Before run(): all accessors that rely on assureFinished() must throw.
    EXPECT_THROW(baseScore.scores(), std::runtime_error);
    EXPECT_THROW(baseScore.score(0), std::runtime_error);
    EXPECT_THROW(baseScore.score(0, 1), std::runtime_error);
    EXPECT_THROW(baseScore.calculate(/*squared*/ true, /*offset*/ 1.0, /*factor*/ 1.0),
                 std::runtime_error);

    // After run(): scores() must be usable (default EdgeScore run() does no preprocessing).
    baseScore.run();
    EXPECT_NO_THROW({
        const auto &scores = baseScore.scores();
        EXPECT_TRUE(scores.empty());
    });
}

TEST_F(EdgeScoreTest, testDerivedScoresReturnsScoreVector) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.indexEdges();

    TestEdgeScoreHelper score(G);
    score.run();

    const auto &scores = score.scores();
    ASSERT_GE(scores.size(), 2u);
    EXPECT_DOUBLE_EQ(scores[0], 1.0);
    EXPECT_DOUBLE_EQ(scores[1], 2.0);
}

TEST_F(EdgeScoreTest, testDerivedScoreByEdgeIdReturnsScore) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.indexEdges();

    TestEdgeScoreHelper score(G);
    score.run();

    EXPECT_DOUBLE_EQ(score.score(0), 1.0);
    EXPECT_DOUBLE_EQ(score.score(1), 2.0);
}

TEST_F(EdgeScoreTest, testDerivedScoreByEdgeNodesReturnsScore) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.indexEdges();

    TestEdgeScoreHelper score(G);
    score.run();

    EXPECT_DOUBLE_EQ(score.score(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(score.score(1, 2), 2.0);
}

TEST_F(EdgeScoreTest, testCalculateThrowsIfEdgesNotIndexed) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    TestEdgeScoreHelper score(G);
    score.run();

    EXPECT_THROW(score.calculate(/*squared*/ true, /*offset*/ 1.0, /*factor*/ 1.0),
                 std::runtime_error);
}

TEST_F(EdgeScoreTest, testCalculateBuildsWeightedGraphWithLinearWeights) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.indexEdges();

    TestEdgeScoreHelper score(G);
    score.run();

    constexpr edgeweight offset = 1.0;
    constexpr edgeweight factor = 2.0;

    Graph W = score.calculate(/*squared*/ false, offset, factor);
    compareTopology(G, W);
    EXPECT_TRUE(W.isWeighted());

    G.forEdges([&](node u, node v, edgeid eid) {
        const edgeweight expected = offset + factor * static_cast<edgeweight>(eid + 1);
        EXPECT_DOUBLE_EQ(W.weight(u, v), expected);
    });
}

TEST_F(EdgeScoreTest, testCalculateBuildsWeightedGraphWithSquaredWeights) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.indexEdges();

    TestEdgeScoreHelper score(G);
    score.run();

    constexpr edgeweight offset = 0.5;
    constexpr edgeweight factor = 1.5;

    Graph W = score.calculate(/*squared*/ true, offset, factor);
    EXPECT_TRUE(W.isWeighted());
    compareTopology(G, W);

    G.forEdges([&](node u, node v, edgeid eid) {
        const edgeweight s = static_cast<edgeweight>(eid + 1);
        const edgeweight expected = offset + factor * s * s;
        EXPECT_DOUBLE_EQ(W.weight(u, v), expected);
    });
}

TEST_F(EdgeScoreTest, testCalculateOnDirectedWeightedGraphFromScores) {
    Graph G(3, true, true);
    G.addEdge(0, 1, 7.0);
    G.addEdge(1, 2, 11.0);
    G.indexEdges();

    TestEdgeScoreHelper score(G);
    score.run();

    constexpr edgeweight offset = 1.333;
    constexpr edgeweight factor = 10.0;
    Graph W = score.calculate(/*squared=*/false, offset, factor);
    EXPECT_TRUE(W.isWeighted());
    compareTopology(G, W);

    G.forEdges([&](node u, node v, edgeid eid) {
        const edgeweight expected =
            offset + factor * static_cast<edgeweight>(eid + 1); // from scoreData[eid]
        EXPECT_DOUBLE_EQ(W.weight(u, v), expected);
    });
}
} // namespace NetworKit
