/*
 * RIGTest.cpp
 *
 *  Created on: Aug 21, 2026
 *      Author: Mikhail Kirilin
 */

#include <algorithm>
#include <memory>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/GlobalState.hpp>
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
namespace {

using IsomorphismDetails::computeOrdering;
using IsomorphismDetails::prepareRISearch;
using IsomorphismDetails::RIDomains;
using IsomorphismDetails::RIImpl;
using IsomorphismDetails::RIOrdering;
using IsomorphismDetails::RISearchSetup;
using IsomorphismDetails::SearchGraph;
using IsomorphismTest::Case;
using IsomorphismTest::Match;
using Semantics = SubgraphIsomorphism::Semantics;
using testing::Contains;
using testing::SizeIs;
using testing::UnorderedElementsAre;
using testing::UnorderedElementsAreArray;

RISearchSetup setupOf(const Graph &pattern, const Graph &target,
                      std::span<const index> patternNodeLabels,
                      std::span<const index> targetNodeLabels, RI::Variant variant) {
    return prepareRISearch(pattern, target, patternNodeLabels, targetNodeLabels, {}, {}, variant,
                           "RI");
}

RISearchSetup setupOf(const Case &testCase, RI::Variant variant) {
    return prepareRISearch(testCase.pattern, testCase.target, testCase.patternNodeLabels,
                           testCase.targetNodeLabels, testCase.patternEdgeLabels,
                           testCase.targetEdgeLabels, variant, "RI");
}

/// Runs the search on unlabelled graphs with the given domains and returns the matches.
std::vector<Match> matchesOf(const SearchGraph &pattern, const SearchGraph &target,
                             const RIOrdering &ordering, const RIDomains &domains) {
    Aux::SignalHandler handler;
    std::vector<Match> matches;
    RIImpl(pattern, target, {}, {}, ordering, domains, Semantics::MONOMORPHISM, handler,
           [&matches](std::span<const node> match) {
               matches.emplace_back(match.begin(), match.end());
               return true;
           })
        .run();

    return matches;
}

index positionOf(const RIOrdering &ordering, node pu) {
    const auto found = std::ranges::find(ordering.order, pu);
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

/// The variants may order the pattern nodes differently but must find the same matches.
class RIGTest : public testing::TestWithParam<RI::Variant> {
protected:
    void TearDown() override { GlobalState::setReceivedSIGINT(false); }
};

INSTANTIATE_TEST_SUITE_P(Variants, RIGTest, testing::Values(RI::Variant::RI, RI::Variant::RI_DS));

/// Checks of the RI-DS domains, which plain RI does not build.
class RIDSGTest : public testing::Test {};

/// Python passes the variant as a plain integer, so the constructor must reject unknown values.
TEST_P(RIGTest, testUnknownVariantIsRejected) {
    const Graph pattern = IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}});
    const Graph target = IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}, {2, 0}});

    EXPECT_THROW(RI algo(pattern, target, static_cast<RI::Variant>(2)), std::invalid_argument);
    EXPECT_THROW(RI algo(pattern, target, static_cast<RI::Variant>(256)), std::invalid_argument);
    EXPECT_NO_THROW(RI algo(pattern, target, GetParam()));
}

TEST_P(RIGTest, testAgreesWithTheReference) {
    const RI::Variant variant = GetParam();
    const auto make = [variant](const Graph &pattern, const Graph &target, Semantics semantics,
                                count maxMatches) -> std::unique_ptr<SubgraphIsomorphism> {
        return std::make_unique<RI>(pattern, target, variant, semantics, maxMatches);
    };

    IsomorphismTest::expectMatchesReference(make);
    IsomorphismTest::expectRespectsMatchCap(make);
}

TEST_P(RIGTest, testOrderingInvariants) {
    for (const Case &testCase : IsomorphismTest::standardCases()) {
        const RISearchSetup setup = setupOf(testCase, GetParam());
        const SearchGraph &pattern = setup.patternGraph;
        const RIOrdering &ordering = setup.ordering;

        ASSERT_THAT(ordering.order, SizeIs(pattern.numberOfNodes())) << "case: " << testCase.name;
        ASSERT_THAT(ordering.parent, SizeIs(ordering.order.size())) << "case: " << testCase.name;

        std::vector<bool> seen(pattern.upperNodeIdBound(), false);
        for (const node pu : ordering.order) {
            ASSERT_LT(pu, pattern.upperNodeIdBound()) << "case: " << testCase.name;
            EXPECT_TRUE(pattern.hasNode(pu))
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
                EXPECT_TRUE(pattern.hasEdge(pp, pu) || pattern.hasEdge(pu, pp))
                    << "case: " << testCase.name << " - parent is not adjacent";
            }

            // A `none` parent means exactly that no earlier position is adjacent.
            bool anyEarlierAdjacent = false;
            for (index j = 0; j < i && !anyEarlierAdjacent; ++j) {
                const node earlier = ordering.order[j];
                anyEarlierAdjacent = pattern.hasEdge(earlier, pu) || pattern.hasEdge(pu, earlier);
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
    // Under RI-DS, the target shapes the order through the domains.
    const auto orderingOf = [variant](const Graph &pattern) {
        return setupOf(pattern, pattern, {}, {}, variant).ordering;
    };

    // The middle node of a path has the maximum degree and goes first.
    const RIOrdering path = orderingOf(IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}}));
    EXPECT_EQ(path.order, (std::vector<node>{1, 0, 2}));
    EXPECT_EQ(path.parent, (std::vector<index>{none, 0, 0}));

    // Position 2 starts a new component.
    const RIOrdering split = orderingOf(IsomorphismTest::graphOf(4, {{0, 1}, {2, 3}}));
    EXPECT_EQ(split.order, (std::vector<node>{0, 1, 2, 3}));
    EXPECT_EQ(split.parent, (std::vector<index>{none, 0, none, 2}));

    const RIOrdering worked = orderingOf(IsomorphismTest::graphOf(
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

    const RISearchSetup tied =
        setupOf(tiedPattern, tiedTarget, tiedPatternLabels, tiedTargetLabels, variant);

    if (variant == RI::Variant::RI) {
        EXPECT_EQ(tied.ordering.order, (std::vector<node>{0, 1}))
            << "with no domains a full tie goes to the smallest node id";
        return;
    }

    ASSERT_THAT(tied.domains.ofPatternNode[0], SizeIs(3));
    ASSERT_THAT(tied.domains.ofPatternNode[1], SizeIs(2));
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

    ASSERT_LT(setupOf(pattern, karate, {}, {}, RI::Variant::RI_DS).domains.ofPatternNode[0].size(),
              karate.numberOfNodes());

    RI plain(pattern, karate, RI::Variant::RI, Semantics::MONOMORPHISM, 0);
    plain.run();

    RI withDomains(pattern, karate, RI::Variant::RI_DS, Semantics::MONOMORPHISM, 0);
    withDomains.run();

    EXPECT_THAT(withDomains.getMatches(), UnorderedElementsAreArray(plain.getMatches()))
        << "RI-DS is pure pruning and must not change the match set";
    EXPECT_THAT(plain.getMatches(), SizeIs(expected));
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
    // Class-0 nodes 1..9 have degree but no class-1 neighbor.
    for (node u = 1; u + 1 < perClass; ++u)
        target.addEdge(u, u + 1);
    target.addEdge(0, perClass);

    std::vector<index> targetNodeLabels(2 * perClass);
    for (node v = 0; v < 2 * perClass; ++v)
        targetNodeLabels[v] = v < perClass ? 0 : 1;

    // Pattern node 1 carries class 0, so the order puts it second, at the position with a parent.
    const Graph pattern = IsomorphismTest::graphOf(2, {{0, 1}});
    const std::vector<index> patternNodeLabels{1, 0};

    const RISearchSetup setup =
        setupOf(pattern, target, patternNodeLabels, targetNodeLabels, RI::Variant::RI_DS);
    ASSERT_EQ(setup.ordering.parent, (std::vector<index>{none, 0}));
    ASSERT_TRUE(setup.domains.earnsItsKeep[1]) << "RI-DS no longer intersects the candidates of "
                                                  "pattern node 1, so this case tests nothing";

    const std::vector<Match> expected = IsomorphismTest::referenceMatches(
        pattern, target, Semantics::MONOMORPHISM, patternNodeLabels, targetNodeLabels);
    ASSERT_FALSE(expected.empty()) << "a case with no matches would not exercise the search";

    RI algo(pattern, target, GetParam(), Semantics::MONOMORPHISM, 0);
    algo.setNodeLabels(patternNodeLabels, targetNodeLabels);
    algo.run();

    EXPECT_THAT(algo.getMatches(), UnorderedElementsAreArray(expected));
}

/// The domains are set by hand, since a computed domain keeps every neighbor that could host the
/// pattern node. The slice of target node 5 has a node below and a node above the domain.
TEST_F(RIDSGTest, testCandidatesAreIntersectedWithTheDomain) {
    const SearchGraph pattern(IsomorphismTest::graphOf(2, {{0, 1}}), /* buildMatrix = */ true);
    const SearchGraph target(IsomorphismTest::graphOf(10, {{5, 1}, {5, 3}, {5, 7}, {5, 9}}),
                             /* buildMatrix = */ false);

    RIDomains domains;
    domains.ofPatternNode = {{5}, {3, 7}};
    domains.earnsItsKeep = {false, true};
    const RIOrdering ordering = computeOrdering(pattern, domains);
    ASSERT_EQ(ordering.parent, (std::vector<index>{none, 0}));

    EXPECT_THAT(matchesOf(pattern, target, ordering, domains),
                UnorderedElementsAre(Match{5, 3}, Match{5, 7}));
}

/// The directed counterpart. Pattern node 1 points to its parent, so its candidates come from the
/// in-slice of target node 5. Domain node 4 is only an out-neighbor of node 5 and must not appear.
TEST_F(RIDSGTest, testInSliceCandidatesAreIntersectedWithTheDomain) {
    const SearchGraph pattern(IsomorphismTest::graphOf(2, {{1, 0}}, /* directed = */ true),
                              /* buildMatrix = */ true);
    const SearchGraph target(IsomorphismTest::graphOf(10, {{1, 5}, {3, 5}, {7, 5}, {9, 5}, {5, 4}},
                                                      /* directed = */ true),
                             /* buildMatrix = */ false);

    RIDomains domains;
    domains.ofPatternNode = {{5}, {3, 4, 7}};
    domains.earnsItsKeep = {false, true};
    const RIOrdering ordering = computeOrdering(pattern, domains);
    ASSERT_EQ(ordering.order, (std::vector<node>{0, 1}));
    ASSERT_EQ(ordering.parent, (std::vector<index>{none, 0}));

    EXPECT_THAT(matchesOf(pattern, target, ordering, domains),
                UnorderedElementsAre(Match{5, 3}, Match{5, 7}));
}

/**
 * A singleton domain opens the matching order, and forward checking removes its target node from
 * all other domains. Pattern node 2 is isolated, so the refinement cannot relate it to nodes 0 and
 * 1. Only forward checking can then remove target node 3 from their domains.
 */
TEST_P(RIGTest, testSingletonDomainsAreRemovedAndComeFirst) {
    const auto expectSingletonsFirst = [](const RISearchSetup &prep, const char *name) {
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

    for (const Case &testCase : IsomorphismTest::standardCases())
        expectSingletonsFirst(setupOf(testCase, GetParam()), testCase.name.c_str());

    constexpr node unique = 3;
    constexpr index rareLabel = 7;

    const Graph target = IsomorphismTest::graphOf(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});
    std::vector<index> targetNodeLabels(6, 0);
    targetNodeLabels[unique] = rareLabel;

    const Graph pattern = IsomorphismTest::graphOf(3, {{0, 1}});
    const std::vector<index> patternNodeLabels{none, none, rareLabel};

    const RISearchSetup prep =
        setupOf(pattern, target, patternNodeLabels, targetNodeLabels, GetParam());
    expectSingletonsFirst(prep, "engineered-singleton");

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
    EXPECT_THAT(domains[1], Contains(2u))
        << "target node 2 is a neighbor of target node 3, which is why the sweep kept 3 - without "
           "this the removal above could be the sweep's rather than forward checking's";

    EXPECT_EQ(prep.ordering.order.front(), 2u)
        << "the singleton-domain node must open the order under RI-DS";
}

/**
 * Removing a singleton's target node can create a new singleton. Degrees give the four star
 * centers nested domains of sizes 1 to 4, so all four end as singletons only if the removals
 * chain.
 */
TEST_F(RIDSGTest, testForwardCheckingReachesAFixpoint) {
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

    // Centers 0..3 have 4, 3, 2 and 1 leaves. The single refinement pass visits the centers before
    // it refines any leaf domain.
    const Graph pattern = IsomorphismTest::graphOf(
        14, {{0, 4}, {0, 5}, {0, 6}, {0, 7}, {1, 8}, {1, 9}, {1, 10}, {2, 11}, {2, 12}, {3, 13}});
    std::vector<index> patternNodeLabels(14, leafLabel);
    for (node u = 0; u < 4; ++u)
        patternNodeLabels[u] = 0;

    const RIDomains domains =
        setupOf(pattern, target, patternNodeLabels, targetNodeLabels, RI::Variant::RI_DS).domains;

    EXPECT_FALSE(domains.anyEmpty);
    for (node center = 0; center < 4; ++center)
        EXPECT_EQ(domains.ofPatternNode[center], (std::vector<node>{center}))
            << "center " << center << " should have been narrowed to one target node by the "
            << "cascade, not left at " << domains.ofPatternNode[center].size() << " entries";
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

    const RIDomains domains =
        setupOf(pattern, target, patternNodeLabels, targetNodeLabels, RI::Variant::RI_DS).domains;
    EXPECT_TRUE(domains.anyEmpty) << "one target node cannot host two pattern nodes at once";

    RI algo(pattern, target, GetParam(), Semantics::MONOMORPHISM, 0);
    algo.setNodeLabels(patternNodeLabels, targetNodeLabels);
    algo.run();
    EXPECT_EQ(algo.numberOfMatches(), 0u);
}

/// The search checks for an interrupt before every report, so no match follows the interrupt.
TEST_P(RIGTest, testInterruptStopsTheSearch) {
    const Graph pattern = IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}});
    const Graph cycle =
        IsomorphismTest::graphOf(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}});

    RI algo(pattern, cycle, GetParam(), Semantics::MONOMORPHISM, 0);
    count delivered = 0;
    algo.setCallback([&](std::span<const node>) {
        if (++delivered == 5)
            GlobalState::setReceivedSIGINT(true);
    });

    EXPECT_THROW(algo.run(), Aux::SignalHandler::InterruptException);
    GlobalState::setReceivedSIGINT(false);

    EXPECT_EQ(delivered, 5u);
    EXPECT_FALSE(algo.hasFinished());
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

                    const std::vector<Match> expected = IsomorphismTest::referenceMatches(
                        pattern, target, semantics, patternNodeLabels, targetNodeLabels);

                    RI algo(pattern, target, variant, semantics, 0);
                    if (nodeLabelled)
                        algo.setNodeLabels(patternNodeLabels, targetNodeLabels);
                    algo.run();

                    EXPECT_THAT(algo.getMatches(), UnorderedElementsAreArray(expected))
                        << "directed: " << directed
                        << ", induced: " << (semantics == Semantics::INDUCED)
                        << ", nodeLabelled: " << nodeLabelled << ", trial: " << trial;
                }
            }
        }
    }
}

} // namespace
} // namespace NetworKit
