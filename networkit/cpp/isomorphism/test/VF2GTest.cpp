/*
 * VF2GTest.cpp
 *
 *  Created on: Aug 17, 2026
 *      Author: Alexandra
 */

#include <algorithm>
#include <iterator>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>
#include <networkit/isomorphism/VF2.hpp>

#include "SubgraphIsomorphismTestUtils.hpp"
#include "../SearchGraph.hpp"

namespace NetworKit {

class VF2GTest : public testing::Test {};

TEST_F(VF2GTest, testExceptionWhenBadlyCollapsed) {

    // Target badly collapsed
    Graph pattern = Graph(2, false, false, true);
    Graph target = Graph(2, false, false, true);

    pattern.addEdge(0, 1);

    target.addEdge(0, 1);
    target.addEdge(0, 1);
    target.addEdge(0, 1);

    std::vector<index> patternLabels = {0};
    std::vector<index> targetLabels = {1, 2, 3};

    EXPECT_EQ(pattern.numberOfEdges(), 1);
    EXPECT_EQ(target.numberOfEdges(), 3);

    VF2 vf = VF2(pattern, target);
    vf.setEdgeLabels(patternLabels, targetLabels);
    try {
        vf.run();
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(),
                     "VF2 does not run if target has unequally-labelled collapsed edges.");
    }

    // Pattern badly collapsed
    Graph pattern2 = Graph(2, false, false, true);
    Graph target2 = Graph(3, false, false, true);

    pattern2.addEdge(0, 1);
    pattern2.addEdge(0, 1);
    pattern2.addEdge(0, 1);

    target2.addEdge(0, 1);
    target2.addEdge(1, 2);
    target2.addEdge(2, 0);

    std::vector<index> patternLabels2 = {0, 1, 2};
    std::vector<index> targetLabels2 = {3, 4, 5};

    EXPECT_EQ(pattern2.numberOfEdges(), 3);
    EXPECT_EQ(target2.numberOfEdges(), 3);

    VF2 vf2 = VF2(pattern2, target2);
    vf2.setEdgeLabels(patternLabels2, targetLabels2);
    try {
        vf2.run();
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(),
                     "VF2 does not run if pattern has unequally-labelled collapsed edges.");
    }

    // Target validly collapsed, only equally labelled edges
    Graph pattern3 = Graph(2, false, false, true);
    Graph target3 = Graph(2, false, false, true);

    pattern3.addEdge(0, 1);

    target3.addEdge(0, 1);
    target3.addEdge(0, 1);
    target3.addEdge(0, 1);

    std::vector<index> patternLabels3 = {0};
    std::vector<index> targetLabels3 = {1, 1, 1};

    EXPECT_EQ(pattern3.numberOfEdges(), 1);
    EXPECT_EQ(target3.numberOfEdges(), 3);

    VF2 vf3 = VF2(pattern3, target3);
    vf3.setEdgeLabels(patternLabels3, targetLabels3);
    EXPECT_NO_THROW(vf3.run());

    // Pattern validly collapsed, only eqally labelled edges
    Graph pattern4 = Graph(2, false, false, true);
    Graph target4 = Graph(3, false, false, true);

    pattern4.addEdge(0, 1);
    pattern4.addEdge(0, 1);
    pattern4.addEdge(0, 1);

    target4.addEdge(0, 1);
    target4.addEdge(1, 2);
    target4.addEdge(2, 0);

    std::vector<index> patternLabels4 = {0, 0, 0};
    std::vector<index> targetLabels4 = {1, 2, 3};

    EXPECT_EQ(pattern4.numberOfEdges(), 3);
    EXPECT_EQ(target4.numberOfEdges(), 3);

    VF2 vf4 = VF2(pattern4, target4);
    vf4.setEdgeLabels(patternLabels4, targetLabels4);
    EXPECT_NO_THROW(vf4.run());

    // Target validly collapsed, only unequally labelled self loops
    Graph pattern5 = Graph(2, false, false, true);
    Graph target5 = Graph(2, false, false, true);

    pattern5.addEdge(0, 1);

    target5.addEdge(0, 0);
    target5.addEdge(0, 0);
    target5.addEdge(0, 0);

    std::vector<index> patternLabels5 = {0};
    std::vector<index> targetLabels5 = {1, 1, 1};

    EXPECT_EQ(pattern5.numberOfEdges(), 1);
    EXPECT_EQ(target5.numberOfEdges(), 3);

    VF2 vf5 = VF2(pattern5, target5);
    vf5.setEdgeLabels(patternLabels5, targetLabels5);
    EXPECT_NO_THROW(vf5.run());
}

TEST_F(VF2GTest, testEdgeLabelMatchingSanityCheck) {

    // Sanity checks for edge label matching

    // Test single edge undirected case, once with matching edge labels and once without
    IsomorphismTest::LabelledGraph ptrnUndir = IsomorphismTest::labelledGraphOf(2, {{0, 1, 0}});
    IsomorphismTest::LabelledGraph trgtMatchedUndir =
        IsomorphismTest::labelledGraphOf(2, {{0, 1, 0}});
    IsomorphismTest::LabelledGraph trgtMismatchedUndir =
        IsomorphismTest::labelledGraphOf(2, {{0, 1, 1}});

    VF2 vfMatchedUndir = VF2(ptrnUndir.G, trgtMatchedUndir.G);
    vfMatchedUndir.setEdgeLabels(ptrnUndir.edgeLabels, trgtMatchedUndir.edgeLabels);
    vfMatchedUndir.run();

    VF2 vfMismatchedUndir = VF2(ptrnUndir.G, trgtMismatchedUndir.G);
    vfMismatchedUndir.setEdgeLabels(ptrnUndir.edgeLabels, trgtMismatchedUndir.edgeLabels);
    vfMismatchedUndir.run();

    EXPECT_TRUE(vfMatchedUndir.hasMatch());
    EXPECT_EQ(vfMatchedUndir.numberOfMatches(), 2);
    EXPECT_EQ(vfMatchedUndir.getMatches(), (std::vector<IsomorphismTest::Match>{{0, 1}, {1, 0}}));

    EXPECT_FALSE(vfMismatchedUndir.hasMatch());

    // Test single edge directed case, once with matching edge labels and once without
    IsomorphismTest::LabelledGraph ptrnDir = IsomorphismTest::labelledGraphOf(2, {{0, 1, 0}}, true);
    IsomorphismTest::LabelledGraph trgtMatchedDir =
        IsomorphismTest::labelledGraphOf(2, {{0, 1, 0}}, true);
    IsomorphismTest::LabelledGraph trgtMismatchedDir =
        IsomorphismTest::labelledGraphOf(2, {{0, 1, 1}}, true);

    VF2 vfMatchedDir = VF2(ptrnDir.G, trgtMatchedDir.G);
    vfMatchedDir.setEdgeLabels(ptrnDir.edgeLabels, trgtMatchedDir.edgeLabels);
    vfMatchedDir.run();

    VF2 vfMismatchedDir = VF2(ptrnDir.G, trgtMismatchedDir.G);
    vfMismatchedDir.setEdgeLabels(ptrnDir.edgeLabels, trgtMismatchedDir.edgeLabels);
    vfMismatchedDir.run();

    EXPECT_TRUE(vfMatchedDir.hasMatch());
    EXPECT_EQ(vfMatchedDir.numberOfMatches(), 1);
    EXPECT_EQ(vfMatchedDir.getMatches(), (std::vector<IsomorphismTest::Match>{{0, 1}}));

    EXPECT_FALSE(vfMismatchedDir.hasMatch());

    // Test single edge undirected case with matching edge labels, node labels only allow one way
    // orientation
    IsomorphismTest::LabelledGraph ptrn1 = IsomorphismTest::labelledGraphOf(2, {{0, 1, 0}});
    IsomorphismTest::LabelledGraph trgt1 = IsomorphismTest::labelledGraphOf(2, {{0, 1, 0}});

    VF2 vf1 = VF2(ptrn1.G, trgt1.G);
    vf1.setNodeLabels({1, 2}, {1, 2});
    vf1.setEdgeLabels(ptrn1.edgeLabels, trgt1.edgeLabels);
    vf1.run();

    EXPECT_TRUE(vf1.hasMatch());
    EXPECT_EQ(vf1.numberOfMatches(), 1);
    EXPECT_EQ(vf1.getMatches(), (std::vector<IsomorphismTest::Match>{{0, 1}}));

    // Test edge matching for square pattern with square containing diagonal line target

    IsomorphismTest::LabelledGraph ptrn2 =
        IsomorphismTest::labelledGraphOf(4, {{0, 1, 0}, {1, 2, 42}, {2, 3, 0}, {3, 0, 42}});
    IsomorphismTest::LabelledGraph trgt2 = IsomorphismTest::labelledGraphOf(
        4, {{0, 1, 42}, {1, 2, 0}, {2, 3, 42}, {3, 0, 0}, {0, 2, 42}});

    VF2 vf2INDUCED = VF2(ptrn2.G, trgt2.G, SubgraphIsomorphism::Semantics::INDUCED);
    vf2INDUCED.setEdgeLabels(ptrn2.edgeLabels, trgt2.edgeLabels);
    vf2INDUCED.run();

    EXPECT_FALSE(vf2INDUCED.hasMatch());

    VF2 vf2MONOMORPHISM = VF2(ptrn2.G, trgt2.G, SubgraphIsomorphism::Semantics::MONOMORPHISM);
    vf2MONOMORPHISM.setEdgeLabels(ptrn2.edgeLabels, trgt2.edgeLabels);
    vf2MONOMORPHISM.run();

    EXPECT_TRUE(vf2MONOMORPHISM.hasMatch());
    EXPECT_EQ(vf2MONOMORPHISM.numberOfMatches(), 4);
    EXPECT_EQ(vf2MONOMORPHISM.getMatches(),
              (std::vector<IsomorphismTest::Match>{
                  {0, 3, 2, 1}, {1, 2, 3, 0}, {2, 1, 0, 3}, {3, 0, 1, 2}}));
}

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

TEST_F(VF2GTest, testTrivialCases) {

    // Empty pattern, undirected
    Graph pattern1 = Graph(0);
    Graph target1 = IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}, {2, 0}});

    VF2 vf1 = VF2(pattern1, target1);
    vf1.run();

    EXPECT_TRUE(vf1.hasMatch());
    EXPECT_EQ(vf1.numberOfMatches(), 1);
    EXPECT_EQ(vf1.getMatches(), std::vector<IsomorphismTest::Match>{{}});

    // Empty pattern, directed
    Graph pattern2 = Graph(0, false, true);
    Graph target2 = IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}, {2, 0}}, true);

    VF2 vf2 = VF2(pattern2, target2);
    vf2.run();

    EXPECT_TRUE(vf2.hasMatch());
    EXPECT_EQ(vf2.numberOfMatches(), 1);
    EXPECT_EQ(vf2.getMatches(), std::vector<IsomorphismTest::Match>{{}});

    // Pattern with more nodes than target, undirected
    Graph pattern3 = Graph(4);
    Graph target3 = Graph(2);

    VF2 vf3 = VF2(pattern3, target3);
    vf3.run();

    EXPECT_FALSE(vf3.hasMatch());
    EXPECT_EQ(vf3.numberOfMatches(), 0);
    EXPECT_EQ(vf3.getMatches(), std::vector<IsomorphismTest::Match>{});

    // Pattern with more nodes than target, directed
    Graph pattern4 = Graph(4, false, true);
    Graph target4 = Graph(2, false, true);

    VF2 vf4 = VF2(pattern4, target4);
    vf4.run();

    EXPECT_FALSE(vf4.hasMatch());
    EXPECT_EQ(vf4.numberOfMatches(), 0);
    EXPECT_EQ(vf4.getMatches(), std::vector<IsomorphismTest::Match>{});

    // Pattern with higher maximum degree than target, undirected
    Graph pattern5 = IsomorphismTest::graphOf(4, {{0, 1}, {0, 2}, {0, 3}});
    Graph target5 = IsomorphismTest::graphOf(5, {{1, 2}});

    VF2 vf5 = VF2(pattern5, target5);
    vf5.run();

    EXPECT_FALSE(vf5.hasMatch());
    EXPECT_EQ(vf5.numberOfMatches(), 0);
    EXPECT_EQ(vf5.getMatches(), std::vector<IsomorphismTest::Match>{});

    // Pattern with higher maximum in-degree than target, directed
    Graph pattern6 = IsomorphismTest::graphOf(4, {{1, 0}, {2, 0}, {3, 0}}, true);
    Graph target6 = IsomorphismTest::graphOf(5, {{1, 2}}, true);

    VF2 vf6 = VF2(pattern6, target6);
    vf6.run();

    EXPECT_FALSE(vf6.hasMatch());
    EXPECT_EQ(vf6.numberOfMatches(), 0);
    EXPECT_EQ(vf6.getMatches(), std::vector<IsomorphismTest::Match>{});

    // Pattern with higher maximum out-degree than target, directed
    Graph pattern7 = IsomorphismTest::graphOf(4, {{0, 1}, {0, 2}, {0, 3}}, true);
    Graph target7 = IsomorphismTest::graphOf(5, {{1, 2}}, true);

    VF2 vf7 = VF2(pattern7, target7);
    vf7.run();

    EXPECT_FALSE(vf7.hasMatch());
    EXPECT_EQ(vf7.numberOfMatches(), 0);
    EXPECT_EQ(vf7.getMatches(), std::vector<IsomorphismTest::Match>{});
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
