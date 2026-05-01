#include <gtest/gtest.h>

#include <networkit/edgescores/SimRankScore.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class SimRankScoreGTest : public ::testing::Test {};

TEST_F(SimRankScoreGTest, testConstructorAcceptsValidParameters) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    EXPECT_NO_THROW({
        SimRankScore score(G, 0.9, 100, 1e-4);
    });
}

TEST_F(SimRankScoreGTest, testConstructorAcceptsDefaultParameters) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    EXPECT_NO_THROW({
        SimRankScore score(G);
    });
}

TEST_F(SimRankScoreGTest, testConstructorDoesNotMarkAlgorithmAsRun) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    SimRankScore score(G);

    EXPECT_THROW(score.scores(), std::runtime_error);
    EXPECT_THROW(score.score(0), std::runtime_error);
    EXPECT_THROW(score.score(0, 1), std::runtime_error);
}

TEST_F(SimRankScoreGTest, testConstructorRejectsNegativeDamping) {
    Graph G(3, false, false);

    EXPECT_THROW(SimRankScore score(G, -0.1, 100, 1e-4), std::invalid_argument);
}

TEST_F(SimRankScoreGTest, testConstructorRejectsDampingGreaterThanOne) {
    Graph G(3, false, false);

    EXPECT_THROW(SimRankScore score(G, 1.1, 100, 1e-4), std::invalid_argument);
}

TEST_F(SimRankScoreGTest, testConstructorRejectsZeroMaxIterations) {
    Graph G(3, false, false);

    EXPECT_THROW(SimRankScore score(G, 0.9, 0, 1e-4), std::invalid_argument);
}

TEST_F(SimRankScoreGTest, testConstructorRejectsNegativeTolerance) {
    Graph G(3, false, false);

    EXPECT_THROW(SimRankScore score(G, 0.9, 100, -1e-4), std::invalid_argument);
}

} // namespace NetworKit