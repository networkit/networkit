/*
 * VF2GTest.cpp
 *
 *  Created on: Aug 17, 2026
 *      Author: Alexandra
 */

#include <memory>

#include <gtest/gtest.h>

#include <networkit/GlobalState.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>
#include <networkit/isomorphism/VF2.hpp>

#include "SubgraphIsomorphismTestUtils.hpp"

namespace NetworKit {

class VF2GTest : public testing::Test {};

TEST_F(VF2GTest, testInitialTriangleCount) {

    // Count triangles in karate graph as initial test
    const Graph pattern = IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}, {2, 0}});

    METISGraphReader reader;
    Graph karate = reader.read("input/karate.graph");

    VF2 vf(pattern, karate);
    vf.run();
    EXPECT_TRUE(vf.hasMatch());
    EXPECT_EQ(vf.numberOfMatches(), 270u);
}

TEST_F(VF2GTest, testAgreesWithTheReference) {

    const auto make = [](const Graph &pattern, const Graph &target,
                         SubgraphIsomorphism::Semantics semantics,
                         count maxMatches) -> std::unique_ptr<SubgraphIsomorphism> {
        return std::make_unique<VF2>(pattern, target, semantics, maxMatches);
    };

    IsomorphismTest::expectMatchesReference(make);
    IsomorphismTest::expectRespectsMatchCap(make);
}

TEST_F(VF2GTest, testInterruptStopsTheSearch) {

    // VF2 polls for CTRL+C once per 1024 candidates, and this search tries about 25000.
    const Graph pattern = IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}});
    const Graph target = IsomorphismTest::completeGraph(30);
    const count total = 30 * 29 * 28;

    VF2 algo(pattern, target, SubgraphIsomorphism::Semantics::MONOMORPHISM);
    count delivered = 0;
    algo.setCallback([&](const SubgraphIsomorphism::Match &) {
        if (++delivered == 1)
            GlobalState::setReceivedSIGINT(true);
    });

    // EXPECT_THROW never ends the test, so the reset always runs.
    EXPECT_THROW(algo.run(), Aux::SignalHandler::InterruptException);
    GlobalState::setReceivedSIGINT(false);

    EXPECT_FALSE(algo.hasFinished());
    EXPECT_LT(delivered, total) << "the search ignored the interrupt";
}

} // namespace NetworKit
