/*
 * RIGTest.cpp
 *
 *  Created on: Aug 21, 2026
 *      Author: Mikhail Kirilin
 */

#include <algorithm>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/edgescores/ChibaNishizekiTriangleEdgeScore.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/isomorphism/RI.hpp>
#include <networkit/isomorphism/SubgraphIsomorphism.hpp>

#include "SubgraphIsomorphismTestUtils.hpp"
#include "../RIImpl.hpp"
#include "../SearchGraph.hpp"

namespace NetworKit {

using IsomorphismDetails::RIImpl;
using IsomorphismDetails::SearchGraph;
using IsomorphismTest::Case;
using IsomorphismTest::Match;
using Semantics = SubgraphIsomorphism::Semantics;

namespace {

struct Snapshot {
    SearchGraph pattern;
    SearchGraph target;

    explicit Snapshot(const Case &testCase)
        : pattern(testCase.pattern, /* buildMatrix = */ true, testCase.patternEdgeLabels),
          target(testCase.target, /* buildMatrix = */ false, testCase.targetEdgeLabels) {}

    bool refused() const {
        return pattern.collapsedLabelledEdges() || target.collapsedLabelledEdges();
    }
};

/// Computes the domains before the ordering, as both drivers do.
struct Preprocessed {
    RIImpl::Domains domains;
    RIImpl::Ordering ordering;

    Preprocessed(const SearchGraph &pattern, const SearchGraph &target,
                 const std::vector<index> &patternNodeLabels,
                 const std::vector<index> &targetNodeLabels, RI::Variant variant)
        : domains(RIImpl::computeDomains(pattern, target, patternNodeLabels, targetNodeLabels,
                                         variant)),
          ordering(RIImpl::computeOrdering(pattern, domains)) {}
};

index positionOf(const RIImpl::Ordering &ordering, node pu) {
    const auto found = std::find(ordering.order.begin(), ordering.order.end(), pu);
    return found == ordering.order.end() ? none
                                         : static_cast<index>(found - ordering.order.begin());
}

std::vector<index> randomNodeLabels(count upperNodeIdBound) {
    std::vector<index> labels;
    labels.reserve(upperNodeIdBound);
    for (count i = 0; i < upperNodeIdBound; ++i) {
        const index drawn = static_cast<index>(Aux::Random::integer(0, 2));
        labels.push_back(drawn == 2 ? none : drawn);
    }
    return labels;
}

} // namespace

/// The variants may order the pattern nodes differently but must find the same matches.
class RIGTest : public testing::TestWithParam<RI::Variant> {};

INSTANTIATE_TEST_SUITE_P(Variants, RIGTest, testing::Values(RI::Variant::RI, RI::Variant::RI_DS));

TEST_P(RIGTest, testAgreesWithTheReference) {
    const RI::Variant variant = GetParam();
    const auto make = [variant](const Graph &pattern, const Graph &target, Semantics semantics,
                                count maxMatches) {
        return std::unique_ptr<SubgraphIsomorphism>(
            new RI(pattern, target, variant, semantics, maxMatches));
    };

    IsomorphismTest::expectMatchesReference(make);
    IsomorphismTest::expectRespectsMatchCap(make);
    IsomorphismTest::expectCallbackFormsAgree(make);
}

TEST_P(RIGTest, testOrderingInvariants) {
    for (const Case &testCase : IsomorphismTest::standardCases()) {
        const Snapshot snapshot(testCase);
        const RIImpl::Ordering ordering =
            Preprocessed(snapshot.pattern, snapshot.target, testCase.patternNodeLabels,
                         testCase.targetNodeLabels, GetParam())
                .ordering;

        ASSERT_EQ(ordering.order.size(), snapshot.pattern.numberOfNodes())
            << "case: " << testCase.name;
        ASSERT_EQ(ordering.parent.size(), ordering.order.size()) << "case: " << testCase.name;

        std::vector<bool> seen(snapshot.pattern.upperNodeIdBound(), false);
        for (const node pu : ordering.order) {
            ASSERT_LT(pu, snapshot.pattern.upperNodeIdBound()) << "case: " << testCase.name;
            EXPECT_TRUE(snapshot.pattern.hasNode(pu))
                << "case: " << testCase.name << " - ordered an id that is not a node";
            EXPECT_FALSE(seen[pu]) << "case: " << testCase.name << " - node ordered twice";
            seen[pu] = true;
        }

        if (!ordering.parent.empty()) {
            EXPECT_EQ(ordering.parent[0], none) << "case: " << testCase.name;
        }

        for (index i = 0; i < ordering.order.size(); ++i) {
            const node pu = ordering.order[i];

            const index parentPos = ordering.parent[i];
            if (parentPos != none) {
                ASSERT_LT(parentPos, i) << "case: " << testCase.name;
                const node pp = ordering.order[parentPos];
                EXPECT_TRUE(snapshot.pattern.hasEdge(pp, pu) || snapshot.pattern.hasEdge(pu, pp))
                    << "case: " << testCase.name << " - parent is not adjacent";
            }

            // A `none` parent means exactly that no earlier position is adjacent.
            bool anyEarlierAdjacent = false;
            for (index j = 0; j < i && !anyEarlierAdjacent; ++j) {
                const node earlier = ordering.order[j];
                anyEarlierAdjacent =
                    snapshot.pattern.hasEdge(earlier, pu) || snapshot.pattern.hasEdge(pu, earlier);
            }
            EXPECT_EQ(parentPos == none, !anyEarlierAdjacent)
                << "case: " << testCase.name << " at position " << i;
        }
    }
}

/**
 * Orders small enough to trace by hand. In the eight-node graph, nodes 5 and 0 tie on the first and
 * third term of the score, so only the two-hop term can order 5 first.
 */
TEST_P(RIGTest, testOrderingHandTraced) {
    const RI::Variant variant = GetParam();
    const auto orderingOf = [variant](const Graph &pattern) {
        const SearchGraph patternGraph(pattern, /* buildMatrix = */ true);
        // Under RI-DS, the target shapes the order through the domains.
        const SearchGraph targetGraph(pattern, /* buildMatrix = */ false);
        return Preprocessed(patternGraph, targetGraph, {}, {}, variant).ordering;
    };

    // The middle node of a path has the maximum degree and goes first.
    const RIImpl::Ordering path = orderingOf(IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}}));
    EXPECT_EQ(path.order, (std::vector<node>{1, 0, 2}));
    EXPECT_EQ(path.parent, (std::vector<index>{none, 0, 0}));

    // Position 2 starts a new component.
    const RIImpl::Ordering split = orderingOf(IsomorphismTest::graphOf(4, {{0, 1}, {2, 3}}));
    EXPECT_EQ(split.order, (std::vector<node>{0, 1, 2, 3}));
    EXPECT_EQ(split.parent, (std::vector<index>{none, 0, none, 2}));

    const RIImpl::Ordering worked = orderingOf(IsomorphismTest::graphOf(
        8, {{4, 1}, {4, 5}, {4, 0}, {4, 3}, {1, 2}, {1, 6}, {5, 2}, {5, 7}, {0, 7}, {3, 6}}));
    EXPECT_EQ(worked.order[0], 4u) << "the first node must be the unique maximum-degree node";
    EXPECT_LT(positionOf(worked, 5), positionOf(worked, 0))
        << "node 5 beats node 0 only on the two-hop term";

    // Two isolated pattern nodes tie on the whole triple. The rarer label gives node 1 the smaller
    // domain without making it a singleton.
    const Graph tiedTarget(5);
    const std::vector<index> tiedTargetLabels{1, 1, 1, 2, 2};

    const Graph tiedPattern(2);
    const std::vector<index> tiedPatternLabels{1, 2};

    const SearchGraph tiedPatternGraph(tiedPattern, /* buildMatrix = */ true);
    const SearchGraph tiedTargetGraph(tiedTarget, /* buildMatrix = */ false);
    const Preprocessed tied(tiedPatternGraph, tiedTargetGraph, tiedPatternLabels, tiedTargetLabels,
                            variant);

    if (variant == RI::Variant::RI) {
        EXPECT_EQ(tied.ordering.order, (std::vector<node>{0, 1}))
            << "with no domains a full tie goes to the smallest node id";
        return;
    }

    ASSERT_EQ(tied.domains.ofPatternNode[0].size(), 3u);
    ASSERT_EQ(tied.domains.ofPatternNode[1].size(), 2u);
    EXPECT_EQ(tied.ordering.order, (std::vector<node>{1, 0}))
        << "the tie on all three counts must go to the more constrained node";
}

/// Unlike on the corpus, the domains prune on karate. ChibaNishizekiTriangleEdgeScore provides an
/// independent expected count.
TEST_P(RIGTest, testVariantsAgreeOnKarate) {
    METISGraphReader reader;
    Graph karate = reader.read("input/karate.graph");
    karate.indexEdges();

    ChibaNishizekiTriangleEdgeScore triangleScore(karate);
    triangleScore.run();
    const std::vector<count> perEdge = triangleScore.scores();
    const count edgeSum = std::accumulate(perEdge.begin(), perEdge.end(), count{0});
    ASSERT_GT(edgeSum, 0u) << "karate should contain triangles";

    // Each triangle is counted on three edges and found once per automorphism, so the number of
    // matches is 6 * (edgeSum / 3) = 2 * edgeSum.
    const count expected = 2 * edgeSum;

    const Graph pattern = IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}, {2, 0}});

    count couldHostATriangleNode = 0;
    karate.forNodes([&](node v) {
        if (karate.degree(v) >= 2)
            ++couldHostATriangleNode;
    });
    ASSERT_LT(couldHostATriangleNode, karate.numberOfNodes());

    RI plain(pattern, karate, RI::Variant::RI, Semantics::MONOMORPHISM, 0);
    plain.run();
    std::vector<Match> plainMatches = plain.getMatches();
    IsomorphismTest::sortMatches(plainMatches);

    RI withDomains(pattern, karate, RI::Variant::RI_DS, Semantics::MONOMORPHISM, 0);
    withDomains.run();
    std::vector<Match> domainMatches = withDomains.getMatches();
    IsomorphismTest::sortMatches(domainMatches);

    EXPECT_EQ(domainMatches, plainMatches) << "RI-Ds is pure pruning and must not change the "
                                              "match set";
    EXPECT_EQ(plainMatches.size(), expected);
    EXPECT_EQ(expected, 270u);
}

/**
 * RI-DS intersects a candidate slice with a domain only if the refinement removed most of the
 * domain, which never happens on the corpus. Here, a single edge joins two classes of ten nodes,
 * so the refinement keeps one of the class-0 nodes.
 */
TEST_P(RIGTest, testSelectiveDomainsDoNotChangeTheMatchSet) {
    constexpr count perClass = 10;

    Graph target(2 * perClass);
    // Class-0 nodes 1..9 have degree but no class-1 neighbour.
    for (node u = 1; u + 1 < perClass; ++u)
        target.addEdge(u, u + 1);
    target.addEdge(0, perClass);

    std::vector<index> targetNodeLabels(2 * perClass);
    for (node v = 0; v < 2 * perClass; ++v)
        targetNodeLabels[v] = v < perClass ? 0 : 1;

    // Pattern node 1 carries class 0, so the order puts it second, at the position with a parent.
    const Graph pattern = IsomorphismTest::graphOf(2, {{0, 1}});
    const std::vector<index> patternNodeLabels{1, 0};

    // At most a fifth of the candidates may survive the refinement, or RI-DS stops intersecting.
    count couldHost = 0;
    count survivesSweep = 0;
    for (node v = 0; v < perClass; ++v) {
        if (target.degree(v) == 0)
            continue;
        ++couldHost;
        target.forNeighborsOf(v, [&](node w) {
            if (targetNodeLabels[w] == 1 && target.degree(w) != 0)
                ++survivesSweep;
        });
    }
    ASSERT_GT(couldHost, 0u);
    ASSERT_LE(survivesSweep * 5, couldHost) << "the sweep is no longer selective enough to make "
                                               "RI-Ds intersect, so this case tests nothing";

    std::vector<Match> expected = IsomorphismTest::referenceMatches(
        pattern, target, Semantics::MONOMORPHISM, patternNodeLabels, targetNodeLabels);
    IsomorphismTest::sortMatches(expected);
    ASSERT_FALSE(expected.empty()) << "a case with no matches would not exercise the search";

    RI algo(pattern, target, GetParam(), Semantics::MONOMORPHISM, 0);
    algo.setNodeLabels(patternNodeLabels, targetNodeLabels);
    algo.run();

    std::vector<Match> actual = algo.getMatches();
    IsomorphismTest::sortMatches(actual);
    EXPECT_EQ(actual, expected);
}

/**
 * A singleton domain opens the matching order, and forward checking removes its target node from
 * all other domains. Pattern node 2 is isolated, so the refinement cannot relate it to nodes 0 and
 * 1. Only forward checking can then remove target node 3 from their domains.
 */
TEST_P(RIGTest, testSingletonDomainsAreRemovedAndComeFirst) {
    const auto expectSingletonsFirst = [](const SearchGraph &pattern, const SearchGraph &target,
                                          const std::vector<index> &patternNodeLabels,
                                          const std::vector<index> &targetNodeLabels,
                                          RI::Variant variant, const char *name) {
        const Preprocessed prep(pattern, target, patternNodeLabels, targetNodeLabels, variant);
        if (prep.domains.ofPatternNode.empty())
            return;

        bool seenLargerDomain = false;
        for (const node pu : prep.ordering.order) {
            if (prep.domains.ofPatternNode[pu].size() == 1)
                EXPECT_FALSE(seenLargerDomain)
                    << "case: " << name << " - pattern node " << pu
                    << " has a singleton domain but is ordered after one that does not";
            else
                seenLargerDomain = true;
        }
    };

    for (const Case &testCase : IsomorphismTest::standardCases()) {
        const Snapshot snapshot(testCase);
        expectSingletonsFirst(snapshot.pattern, snapshot.target, testCase.patternNodeLabels,
                              testCase.targetNodeLabels, GetParam(), testCase.name.c_str());
    }

    constexpr node unique = 3;
    constexpr index rareLabel = 7;

    const Graph target = IsomorphismTest::graphOf(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});
    std::vector<index> targetNodeLabels(6, 0);
    targetNodeLabels[unique] = rareLabel;

    const Graph pattern = IsomorphismTest::graphOf(3, {{0, 1}});
    const std::vector<index> patternNodeLabels{none, none, rareLabel};

    const SearchGraph patternGraph(pattern, /* buildMatrix = */ true);
    const SearchGraph targetGraph(target, /* buildMatrix = */ false);

    expectSingletonsFirst(patternGraph, targetGraph, patternNodeLabels, targetNodeLabels,
                          GetParam(), "engineered-singleton");

    const Preprocessed prep(patternGraph, targetGraph, patternNodeLabels, targetNodeLabels,
                            GetParam());

    if (GetParam() == RI::Variant::RI) {
        EXPECT_TRUE(prep.domains.ofPatternNode.empty()) << "plain RI must build no domains at all";
        EXPECT_EQ(prep.ordering.order.front(), 0u)
            << "plain RI has no domains, so it starts at a maximum-degree node as it always did";
        return;
    }

    const std::vector<std::vector<node>> &domains = prep.domains.ofPatternNode;

    EXPECT_FALSE(prep.domains.anyEmpty);
    EXPECT_EQ(domains[2], (std::vector<node>{unique}));
    EXPECT_EQ(domains[0], (std::vector<node>{0, 1, 2, 4, 5}));
    EXPECT_EQ(domains[1], (std::vector<node>{0, 1, 2, 4, 5}));
    EXPECT_TRUE(std::find(domains[1].begin(), domains[1].end(), 2u) != domains[1].end())
        << "target node 2 is a neighbour of target node 3, which is why the sweep kept 3 - without "
           "this the removal above could be the sweep's rather than forward checking's";

    EXPECT_EQ(prep.ordering.order.front(), 2u)
        << "the singleton-domain node must open the order under RI-Ds";
}

/**
 * Removing a singleton's target node can create a new singleton. Degrees give the four star
 * centres nested domains of sizes 1 to 4, so all four end as singletons only if the removals
 * chain.
 */
TEST_P(RIGTest, testForwardCheckingReachesAFixpoint) {
    constexpr index leafLabel = 9;

    // Target nodes 0..3 carry label 0 and have degrees 4, 3, 2 and 1; 4..13 are label-9 pendants.
    const Graph target = IsomorphismTest::graphOf(
        14, {{0, 4}, {0, 11}, {0, 12}, {0, 13}, {1, 5}, {1, 9}, {1, 10}, {2, 6}, {2, 8}, {3, 7}});
    std::vector<index> targetNodeLabels(14, leafLabel);
    for (node v = 0; v < 4; ++v)
        targetNodeLabels[v] = 0;

    for (const count degree : {4u, 3u, 2u, 1u}) {
        count labelZeroNodesOfThatDegree = 0;
        target.forNodes([&](node v) {
            if (targetNodeLabels[v] == 0 && target.degree(v) >= degree)
                ++labelZeroNodesOfThatDegree;
        });
        ASSERT_EQ(labelZeroNodesOfThatDegree, 5 - degree);
    }

    // Centres 0..3 have 4, 3, 2 and 1 leaves. The single refinement pass visits the centres before
    // it refines any leaf domain.
    const Graph pattern = IsomorphismTest::graphOf(
        14, {{0, 4}, {0, 5}, {0, 6}, {0, 7}, {1, 8}, {1, 9}, {1, 10}, {2, 11}, {2, 12}, {3, 13}});
    std::vector<index> patternNodeLabels(14, leafLabel);
    for (node u = 0; u < 4; ++u)
        patternNodeLabels[u] = 0;

    const SearchGraph patternGraph(pattern, /* buildMatrix = */ true);
    const SearchGraph targetGraph(target, /* buildMatrix = */ false);
    const RIImpl::Domains domains = RIImpl::computeDomains(
        patternGraph, targetGraph, patternNodeLabels, targetNodeLabels, GetParam());

    if (GetParam() == RI::Variant::RI) {
        EXPECT_TRUE(domains.ofPatternNode.empty()) << "plain RI must build no domains at all";
        return;
    }

    EXPECT_FALSE(domains.anyEmpty);
    for (node centre = 0; centre < 4; ++centre)
        EXPECT_EQ(domains.ofPatternNode[centre], (std::vector<node>{centre}))
            << "centre " << centre << " should have been narrowed to one target node by the "
            << "cascade, not left at " << domains.ofPatternNode[centre].size() << " entries";
}

/// Two pattern nodes that need the same target node make the instance impossible. Plain RI, which
/// has no domains, confirms that no match exists.
TEST_P(RIGTest, testForwardCheckingRejectsImpossibleInstances) {
    constexpr index rareLabel = 7;

    const Graph target = IsomorphismTest::graphOf(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});
    std::vector<index> targetNodeLabels(6, 0);
    targetNodeLabels[2] = rareLabel;

    // Nodes 0 and 2 lie in different components, so the refinement never relates them.
    const Graph pattern = IsomorphismTest::graphOf(4, {{0, 1}, {2, 3}});
    const std::vector<index> patternNodeLabels{rareLabel, 0, rareLabel, 0};

    const SearchGraph patternGraph(pattern, /* buildMatrix = */ true);
    const SearchGraph targetGraph(target, /* buildMatrix = */ false);
    const RIImpl::Domains domains = RIImpl::computeDomains(
        patternGraph, targetGraph, patternNodeLabels, targetNodeLabels, RI::Variant::RI_DS);
    EXPECT_TRUE(domains.anyEmpty) << "one target node cannot host two pattern nodes at once";

    RI algo(pattern, target, GetParam(), Semantics::MONOMORPHISM, 0);
    algo.setNodeLabels(patternNodeLabels, targetNodeLabels);
    algo.run();
    EXPECT_EQ(algo.numberOfMatches(), 0u);
}

/// Drives RIImpl one level at a time, as the ParallelRI workers do. A difference to run() means
/// that a State lacks something the recursion keeps on its stack.
TEST_P(RIGTest, testExpandAgreesWithRun) {
    const RI::Variant variant = GetParam();

    for (const Case &testCase : IsomorphismTest::standardCases()) {
        const Snapshot snapshot(testCase);
        if (snapshot.refused())
            continue;

        Aux::SignalHandler handler;
        const Preprocessed prep(snapshot.pattern, snapshot.target, testCase.patternNodeLabels,
                                testCase.targetNodeLabels, variant);
        const RIImpl::Ordering &ordering = prep.ordering;

        std::vector<Match> viaRun;
        RIImpl(snapshot.pattern, snapshot.target, testCase.patternNodeLabels,
               testCase.targetNodeLabels, ordering, prep.domains, testCase.semantics, handler,
               [&viaRun](const Match &match) {
                   viaRun.push_back(match);
                   return true;
               })
            .run();

        std::vector<Match> viaExpand;
        RIImpl expander(snapshot.pattern, snapshot.target, testCase.patternNodeLabels,
                        testCase.targetNodeLabels, ordering, prep.domains, testCase.semantics,
                        handler, [&viaExpand](const Match &match) {
                            viaExpand.push_back(match);
                            return true;
                        });

        std::vector<RIImpl::State> pending{expander.rootState()};
        while (!pending.empty()) {
            RIImpl::State state = std::move(pending.back());
            pending.pop_back();

            std::vector<RIImpl::State> children;
            EXPECT_TRUE(expander.expand(state, children)) << "case: " << testCase.name;

            for (RIImpl::State &child : children)
                pending.push_back(std::move(child));
        }

        IsomorphismTest::sortMatches(viaRun);
        IsomorphismTest::sortMatches(viaExpand);
        EXPECT_EQ(viaExpand, viaRun) << "case: " << testCase.name;
    }
}

TEST_P(RIGTest, testMatchesReferenceOnRandomGraphs) {
    const RI::Variant variant = GetParam();

    // The reference enumerates all 6^4 mappings per trial, which keeps the sizes small.
    constexpr count patternNodes = 4;
    constexpr count targetNodes = 6;
    constexpr int trials = 8;

    Aux::Random::setSeed(1701, false);

    for (const bool directed : {false, true}) {
        for (const Semantics semantics : {Semantics::MONOMORPHISM, Semantics::INDUCED}) {
            for (const bool nodeLabelled : {false, true}) {
                for (int trial = 0; trial < trials; ++trial) {
                    const Graph pattern =
                        ErdosRenyiGenerator(patternNodes, 0.5, directed).generate();
                    const Graph target = ErdosRenyiGenerator(targetNodes, 0.4, directed).generate();

                    std::vector<index> patternNodeLabels;
                    std::vector<index> targetNodeLabels;
                    if (nodeLabelled) {
                        patternNodeLabels = randomNodeLabels(pattern.upperNodeIdBound());
                        targetNodeLabels = randomNodeLabels(target.upperNodeIdBound());
                    }

                    std::vector<Match> expected = IsomorphismTest::referenceMatches(
                        pattern, target, semantics, patternNodeLabels, targetNodeLabels);
                    IsomorphismTest::sortMatches(expected);

                    RI algo(pattern, target, variant, semantics, 0);
                    if (nodeLabelled)
                        algo.setNodeLabels(patternNodeLabels, targetNodeLabels);
                    algo.run();

                    std::vector<Match> actual = algo.getMatches();
                    IsomorphismTest::sortMatches(actual);

                    EXPECT_EQ(actual, expected)
                        << "directed: " << directed
                        << ", induced: " << (semantics == Semantics::INDUCED)
                        << ", nodeLabelled: " << nodeLabelled << ", trial: " << trial;
                }
            }
        }
    }
}

} // namespace NetworKit
