/*
 * VF2GTest.cpp
 *
 *  Created on: Aug 17, 2026
 *      Author: Alexandra
 */

#include <memory>

#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>
#include <networkit/isomorphism/VF2.hpp>

#include "SubgraphIsomorphismTestUtils.hpp"

namespace NetworKit {

class VF2GTest : public testing::Test {};

TEST_F(VF2GTest, testInitialTriangleCount) {

    // Count triangles in karate graph as initial test
    Graph pattern = Graph(3);

    pattern.addEdge(0, 1);
    pattern.addEdge(1, 2);
    pattern.addEdge(2, 0);

    METISGraphReader reader;
    Graph karate = reader.read("input/karate.graph");

    VF2 vf = VF2(pattern, karate);
    vf.run();
    EXPECT_TRUE(vf.hasMatch());
    EXPECT_EQ(vf.numberOfMatches(), 270u);
}

TEST_F(VF2GTest, testAgreesWithTheReference) {

    const auto make = [](const Graph &pattern, const Graph &target,
                         SubgraphIsomorphism::Semantics semantics, count maxMatches) {
        return std::unique_ptr<SubgraphIsomorphism>(
            new VF2(pattern, target, semantics, maxMatches));
    };

    IsomorphismTest::expectMatchesReference(make);
    IsomorphismTest::expectRespectsMatchCap(make);
    IsomorphismTest::expectCallbackFormsAgree(make);
}

} // namespace NetworKit
