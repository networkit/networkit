/*
 * SubgraphIsomorphismGTest.cpp
 *
 * Tests the base class SubgraphIsomorphism: the run() protocol, input validation and
 * interruption. Most tests use ReferenceSubgraphIsomorphism, a stand-in that reports a precomputed
 * match set.
 */

#include <memory>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/GlobalState.hpp>
#include <networkit/Globals.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

#include "SubgraphIsomorphismTestUtils.hpp"

namespace NetworKit {

using IsomorphismTest::completeGraph;
using IsomorphismTest::graphOf;
using IsomorphismTest::Match;
using IsomorphismTest::referenceMatches;
using IsomorphismTest::Semantics;

class SubgraphIsomorphismGTest : public testing::Test {
protected:
    /// The interrupt tests set a process-global flag, which a failing ASSERT_* would leave set.
    void TearDown() override { GlobalState::setReceivedSIGINT(false); }
};

TEST_F(SubgraphIsomorphismGTest, testHarnessDrivesAnAlgorithmThroughTheRunProtocol) {

    const auto make = [](const Graph &pattern, const Graph &target, Semantics semantics,
                         count maxMatches) -> std::unique_ptr<SubgraphIsomorphism> {
        return std::make_unique<IsomorphismTest::ReferenceSubgraphIsomorphism>(
            pattern, target, semantics, maxMatches);
    };

    IsomorphismTest::expectMatchesReference(make);
    IsomorphismTest::expectRespectsMatchCap(make);
}

TEST_F(SubgraphIsomorphismGTest, testInterruptLeavesTheAlgorithmUnfinishedButUsable) {

    // The test sets the global flag instead of raising a real SIGINT.
    const Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    const Graph target = completeGraph(6);

    std::vector<Match> expected = referenceMatches(pattern, target, Semantics::MONOMORPHISM);
    IsomorphismTest::sortMatches(expected);

    IsomorphismTest::ReferenceSubgraphIsomorphism algo(pattern, target, Semantics::MONOMORPHISM);

    count delivered = 0;
    algo.setCallback([&](const Match &) {
        if (++delivered == 5)
            GlobalState::setReceivedSIGINT(true);
    });

    EXPECT_THROW(algo.run(), Aux::SignalHandler::InterruptException);
    GlobalState::setReceivedSIGINT(false);

    EXPECT_FALSE(algo.hasFinished()) << "an interrupted run must not count as finished";
    EXPECT_THROW(algo.numberOfMatches(), std::runtime_error);
    EXPECT_THROW(algo.hasMatch(), std::runtime_error);

    EXPECT_GE(delivered, 5u);

    // A second object without a callback, since the recovery check needs the stored matches.
    IsomorphismTest::ReferenceSubgraphIsomorphism recovering(pattern, target,
                                                             Semantics::MONOMORPHISM);

    GlobalState::setReceivedSIGINT(true);
    EXPECT_THROW(recovering.run(), Aux::SignalHandler::InterruptException);
    GlobalState::setReceivedSIGINT(false);
    ASSERT_FALSE(recovering.hasFinished());

    recovering.run();

    ASSERT_TRUE(recovering.hasFinished());
    std::vector<Match> actual = recovering.getMatches();
    IsomorphismTest::sortMatches(actual);
    EXPECT_EQ(actual, expected);
}

TEST_F(SubgraphIsomorphismGTest, testSetEdgeLabelsValidatesItsInput) {

    Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    Graph target = graphOf(4, {{0, 1}, {1, 2}, {2, 3}});

    IsomorphismTest::ReferenceSubgraphIsomorphism unindexed(pattern, target,
                                                            Semantics::MONOMORPHISM);
    ASSERT_FALSE(pattern.hasEdgeIds());
    EXPECT_THROW(unindexed.setEdgeLabels({1, 2}, {1, 2, 3}), std::runtime_error)
        << "edge labels are indexed by edge id, so a graph without ids has no index space";

    // Two empty vectors clear the labels, even without edge ids.
    EXPECT_NO_THROW(unindexed.setEdgeLabels({}, {}));

    pattern.indexEdges();
    EXPECT_THROW(unindexed.setEdgeLabels({1, 2}, {1, 2, 3}), std::runtime_error)
        << "the target needs edge ids as well";

    target.indexEdges();
    IsomorphismTest::ReferenceSubgraphIsomorphism algo(pattern, target, Semantics::MONOMORPHISM);

    EXPECT_THROW(algo.setEdgeLabels({1}, {1, 2, 3}), std::runtime_error) << "pattern too short";
    EXPECT_THROW(algo.setEdgeLabels({1, 2}, {1, 2}), std::runtime_error) << "target too short";
    EXPECT_NO_THROW(algo.setEdgeLabels({1, 2}, {1, 2, 3}));
    EXPECT_NO_THROW(algo.setEdgeLabels({}, {}));
}

} // namespace NetworKit
