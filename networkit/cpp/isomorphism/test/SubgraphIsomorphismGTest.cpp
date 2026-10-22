/*
 * SubgraphIsomorphismGTest.cpp
 *
 * Tests the base class SubgraphIsomorphism: the run() protocol, input validation and
 * interruption. Most tests use ReferenceSubgraphIsomorphism, a stand-in that reports a precomputed
 * match set.
 */

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/GlobalState.hpp>
#include <networkit/Globals.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>
#include <networkit/isomorphism/VF2.hpp>

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

TEST_F(SubgraphIsomorphismGTest, testGraphsMustAgreeOnDirectedness) {

    const Graph undirected = graphOf(2, {{0, 1}});
    const Graph directed = graphOf(2, {{0, 1}}, true);

    EXPECT_THROW(VF2 algo(undirected, directed), std::runtime_error);
    EXPECT_THROW(VF2 algo(directed, undirected), std::runtime_error);
    EXPECT_NO_THROW(VF2 algo(directed, directed));
}

TEST_F(SubgraphIsomorphismGTest, testUnknownSemanticsIsRejected) {

    // Python passes the semantics as a plain integer, so any value can arrive here.
    const Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    const Graph target = graphOf(3, {{0, 1}, {1, 2}, {2, 0}});

    EXPECT_THROW(VF2 algo(pattern, target, static_cast<Semantics>(2)), std::invalid_argument);
    EXPECT_THROW(VF2 algo(pattern, target, static_cast<Semantics>(256)), std::invalid_argument);
    EXPECT_NO_THROW(VF2 algo(pattern, target, Semantics::MONOMORPHISM));
}

TEST_F(SubgraphIsomorphismGTest, testGetMatchesFollowsTheLastRun) {

    const Graph pattern = graphOf(2, {{0, 1}});
    const Graph target = graphOf(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});

    VF2 algo(pattern, target, Semantics::MONOMORPHISM);
    algo.setStoreMatches(false);
    algo.run();
    ASSERT_EQ(algo.numberOfMatches(), 8u);

    // The run only counted, so turning storage on afterwards must not expose an empty result.
    algo.setStoreMatches(true);
    EXPECT_THROW(algo.getMatches(), std::runtime_error);

    algo.run();
    EXPECT_EQ(algo.getMatches().size(), 8u);

    // A callback set after the run does not take away the matches that run stored.
    algo.setCallback([](const Match &) {});
    EXPECT_EQ(algo.getMatches().size(), 8u);
}

TEST_F(SubgraphIsomorphismGTest, testSetNodeLabelsValidatesItsInput) {

    const Graph pattern = graphOf(2, {{0, 1}});
    const Graph target = graphOf(3, {{0, 1}, {1, 2}});
    VF2 algo(pattern, target, Semantics::MONOMORPHISM);

    EXPECT_THROW(algo.setNodeLabels({1}, {1, 1, 1}), std::runtime_error) << "pattern too short";
    EXPECT_THROW(algo.setNodeLabels({1, 1}, {1, 1}), std::runtime_error) << "target too short";
    EXPECT_NO_THROW(algo.setNodeLabels({1, 1}, {1, 1, 1}));
    EXPECT_NO_THROW(algo.setNodeLabels({}, {}));
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

TEST_F(SubgraphIsomorphismGTest, testParallelEdgesWithDisagreeingLabelsAreRefused) {

    // One arc of the snapshot cannot represent the two differently labelled 0-1 edges.
    const IsomorphismTest::LabelledGraph plain =
        IsomorphismTest::labelledGraphOf(4, {{0, 1, 1}, {1, 2, 2}, {2, 3, 1}});
    const IsomorphismTest::LabelledGraph disagreeing =
        IsomorphismTest::labelledGraphOf(4, {{0, 1, 1}, {0, 1, 4}, {1, 2, 2}, {2, 3, 1}});

    for (const std::string side : {"pattern", "target"}) {
        const IsomorphismTest::LabelledGraph &pattern = side == "pattern" ? disagreeing : plain;
        const IsomorphismTest::LabelledGraph &target = side == "pattern" ? plain : disagreeing;

        VF2 algo(pattern.G, target.G, Semantics::MONOMORPHISM);
        algo.setEdgeLabels(pattern.edgeLabels, target.edgeLabels);
        try {
            algo.run();
            ADD_FAILURE() << "disagreeing labels in the " << side << " were accepted";
        } catch (const std::runtime_error &e) {
            EXPECT_NE(std::string(e.what()).find(side), std::string::npos) << e.what();
        }
        EXPECT_FALSE(algo.hasFinished());
    }
}

TEST_F(SubgraphIsomorphismGTest, testRunRechecksGraphInvariants) {

    Graph pattern = graphOf(3, {{0, 1}, {1, 2}});
    const Graph target = graphOf(4, {{0, 1}, {1, 2}, {2, 3}});

    VF2 algo(pattern, target, Semantics::MONOMORPHISM);
    ASSERT_NO_THROW(algo.run());

    pattern.addEdge(1, 1);
    ASSERT_EQ(pattern.numberOfSelfLoops(), 1u);

    EXPECT_THROW(algo.run(), std::runtime_error)
        << "a self-loop added to the pattern after construction must still be rejected";
    EXPECT_FALSE(algo.hasFinished());
}

TEST_F(SubgraphIsomorphismGTest, testRunRechecksNodeLabelSizes) {

    const Graph pattern = graphOf(2, {{0, 1}});
    Graph target = graphOf(3, {{0, 1}, {1, 2}});

    VF2 algo(pattern, target, Semantics::MONOMORPHISM);
    algo.setNodeLabels({1, 1}, {1, 1, 1});
    ASSERT_NO_THROW(algo.run());
    ASSERT_GT(algo.numberOfMatches(), 0u);

    target.addNode();
    ASSERT_GT(target.upperNodeIdBound(), 3u);

    EXPECT_THROW(algo.run(), std::runtime_error)
        << "a node added after setNodeLabels() leaves the label vector short";

    // A rejected run must not hand back the matches of the previous run.
    EXPECT_FALSE(algo.hasFinished());
    EXPECT_THROW(algo.getMatches(), std::runtime_error);

    algo.setNodeLabels({1, 1}, {1, 1, 1, 1});
    EXPECT_NO_THROW(algo.run());
    EXPECT_TRUE(algo.hasFinished());
}

TEST_F(SubgraphIsomorphismGTest, testRunRechecksEdgeLabelSizes) {

    const IsomorphismTest::LabelledGraph pattern = IsomorphismTest::labelledGraphOf(2, {{0, 1, 1}});
    IsomorphismTest::LabelledGraph target =
        IsomorphismTest::labelledGraphOf(3, {{0, 1, 1}, {1, 2, 1}});

    VF2 algo(pattern.G, target.G, Semantics::MONOMORPHISM);
    algo.setEdgeLabels(pattern.edgeLabels, target.edgeLabels);
    ASSERT_NO_THROW(algo.run());

    const index boundBefore = target.G.upperEdgeIdBound();
    target.G.addEdge(0, 2);
    target.G.indexEdges(true);
    ASSERT_GT(target.G.upperEdgeIdBound(), boundBefore);

    EXPECT_THROW(algo.run(), std::runtime_error)
        << "an edge added after setEdgeLabels() leaves the edge label vector short";
    EXPECT_FALSE(algo.hasFinished());

    std::vector<index> grown = target.edgeLabels;
    grown.resize(target.G.upperEdgeIdBound(), 1);
    algo.setEdgeLabels(pattern.edgeLabels, grown);
    EXPECT_NO_THROW(algo.run());
    EXPECT_TRUE(algo.hasFinished());
}

} // namespace NetworKit
