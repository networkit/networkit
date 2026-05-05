#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <networkit/edgescores/SimRankScore.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class SimRankScoreGTest : public ::testing::Test {};

TEST_F(SimRankScoreGTest, testConstructorAcceptsValidParameters) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    EXPECT_NO_THROW({ SimRankScore score(G, 0.9, 100, 1e-4); });
}

TEST_F(SimRankScoreGTest, testConstructorAcceptsDefaultParameters) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    EXPECT_NO_THROW({ SimRankScore score(G); });
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

TEST_F(SimRankScoreGTest, testRunRequiresIndexedEdges) {
    Graph G(2, false, false);
    G.addEdge(0, 1);

    SimRankScore simrank(G);

    EXPECT_THROW(simrank.run(), std::runtime_error);
}

TEST_F(SimRankScoreGTest, testScoresUnavailableBeforeRun) {
    Graph G(2, false, false);
    G.addEdge(0, 1);
    G.indexEdges();

    SimRankScore simrank(G);

    EXPECT_ANY_THROW(simrank.scores());
    EXPECT_ANY_THROW(simrank.score(0));
    EXPECT_ANY_THROW(simrank.score(0, 1));
}

TEST_F(SimRankScoreGTest, testRunProducesOneScorePerIndexedEdge) {
    Graph G(3, false, false, true);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    SimRankScore simrank(G);
    simrank.run();

    EXPECT_THAT(simrank.scores(), testing::SizeIs(G.upperEdgeIdBound()));
}

TEST_F(SimRankScoreGTest, testSingleUndirectedEdgeHasZeroSimRankScore) {
    Graph G(2, false, false);
    G.addEdge(0, 1);
    G.indexEdges();

    SimRankScore simrank(G, 0.9, 100, 1e-12);
    simrank.run();

    EXPECT_NEAR(simrank.score(0, 1), 0.0, 1e-12);
}

TEST_F(SimRankScoreGTest, testPathOfLengthTwoHasZeroEdgeScores) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.indexEdges();

    SimRankScore simrank(G, 0.9, 100, 1e-12);
    simrank.run();

    EXPECT_NEAR(simrank.score(0, 1), 0.0, 1e-12);
    EXPECT_NEAR(simrank.score(1, 2), 0.0, 1e-12);
}

TEST_F(SimRankScoreGTest, testTriangleHasPositiveEqualEdgeScores) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(0, 2);
    G.indexEdges();

    SimRankScore simrank(G, 0.5, 100, 1e-12);
    simrank.run();

    const double expected = 0.2;

    EXPECT_NEAR(simrank.score(0, 1), expected, 1e-10);
    EXPECT_NEAR(simrank.score(1, 2), expected, 1e-10);
    EXPECT_NEAR(simrank.score(0, 2), expected, 1e-10);
}

TEST_F(SimRankScoreGTest, testScoresAreStableAcrossRepeatedRun) {
    Graph G(3, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(0, 2);
    G.indexEdges();

    SimRankScore simrank(G, 0.5, 100, 1e-12);

    simrank.run();
    const auto firstScores = simrank.scores();

    simrank.run();
    const auto secondScores = simrank.scores();

  EXPECT_THAT(firstScores, testing::Pointwise(testing::DoubleNear(1e-12), secondScores));
}

TEST_F(SimRankScoreGTest, testTriangleWithTailHasExpectedEdgeScores) {
    Graph G(4, false, false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(0, 2);
    G.addEdge(2, 3);
    G.indexEdges();

    SimRankScore simrank(G, 0.5, 100, 1e-12);
    simrank.run();

    EXPECT_NEAR(simrank.score(0, 1), 79.0 / 423.0, 1e-10);
    EXPECT_NEAR(simrank.score(1, 2), 65.0 / 423.0, 1e-10);
    EXPECT_NEAR(simrank.score(0, 2), 65.0 / 423.0, 1e-10);
    EXPECT_NEAR(simrank.score(2, 3), 26.0 / 423.0, 1e-10);
}

TEST_F(SimRankScoreGTest, testDirectedGraphUsesInNeighbors) {
    Graph G(3, false, true);

    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(1, 2);

    G.indexEdges();

    SimRankScore simRank(G, 0.8, 1, 0.0);
    simRank.run();

    EXPECT_NEAR(simRank.score(G.edgeId(1, 2)), 0.4, 1e-12);
}

TEST_F(SimRankScoreGTest, testResultsDoNotDependOnNumberOfThreads) {
#ifndef _OPENMP
    GTEST_SKIP() << "OpenMP is not enabled";
#else
    const int previousThreads = omp_get_max_threads();

    Graph G(6, false, false);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(1, 5);
    G.addEdge(0, 4);
    G.indexEdges();

    omp_set_num_threads(1);
    SimRankScore singleThreaded(G, 0.7, 100, 1e-12);
    singleThreaded.run();
    const auto singleThreadedScores = singleThreaded.scores();

    omp_set_num_threads(4);
    SimRankScore multiThreaded(G, 0.7, 100, 1e-12);
    multiThreaded.run();
    const auto multiThreadedScores = multiThreaded.scores();

    omp_set_num_threads(previousThreads);

    EXPECT_EQ(singleThreadedScores.size(), multiThreadedScores.size());

    for (index i = 0; i < singleThreadedScores.size(); ++i) {
        EXPECT_NEAR(singleThreadedScores[i], multiThreadedScores[i], 1e-12);
    }
#endif
}

TEST_F(SimRankScoreGTest, testRepeatedParallelRunsProduceSameScores) {
#ifdef _OPENMP
    const int previousThreads = omp_get_max_threads();
    omp_set_num_threads(4);
#endif

    Graph G(6, false, false);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(3, 4);
    G.addEdge(4, 5);
    G.addEdge(1, 5);
    G.addEdge(0, 4);
    G.indexEdges();

    SimRankScore reference(G, 0.7, 100, 1e-12);
    reference.run();
    const auto referenceScores = reference.scores();

    for (index run = 0; run < 20; ++run) {
        SimRankScore simrank(G, 0.7, 100, 1e-12);
        simrank.run();

        EXPECT_EQ(referenceScores.size(), simrank.scores().size());

        for (index i = 0; i < referenceScores.size(); ++i) {
            EXPECT_NEAR(referenceScores[i], simrank.scores()[i], 1e-12);
        }
    }

#ifdef _OPENMP
    omp_set_num_threads(previousThreads);
#endif
}
} // namespace NetworKit
