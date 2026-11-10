/*
 * RIGTest.cpp
 *
 *  Created on: Aug 21, 2026
 *      Author: Mikhail Kirilin
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
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

using IsomorphismDetails::orderingFrom;
using IsomorphismDetails::prepareRISearch;
using IsomorphismDetails::RIDomains;
using IsomorphismDetails::RIImpl;
using IsomorphismDetails::RIOrdering;
using IsomorphismDetails::RISearchSetup;
using IsomorphismDetails::RIState;
using IsomorphismDetails::SearchGraph;
using IsomorphismTest::Case;
using IsomorphismTest::Match;
using Semantics = SubgraphIsomorphism::Semantics;
using testing::ElementsAre;
using testing::ElementsAreArray;
using testing::HasSubstr;
using testing::IsEmpty;
using testing::Lt;
using testing::Not;
using testing::SizeIs;
using testing::ThrowsMessage;
using testing::UnorderedElementsAre;
using testing::UnorderedElementsAreArray;
using testing::WhenSorted;

/// The default of the RI constructor.
constexpr double defaultDomainFilterThreshold = 0.8;

RISearchSetup setupOf(const Graph &pattern, const Graph &target,
                      std::span<const index> patternNodeLabels,
                      std::span<const index> targetNodeLabels, RI::Variant variant) {
    return prepareRISearch(pattern, target, patternNodeLabels, targetNodeLabels, {}, {}, variant,
                           defaultDomainFilterThreshold, "RI");
}

RISearchSetup setupOf(const Case &testCase, RI::Variant variant,
                      double domainFilterThreshold = defaultDomainFilterThreshold) {
    return prepareRISearch(testCase.pattern, testCase.target, testCase.patternNodeLabels,
                           testCase.targetNodeLabels, testCase.patternEdgeLabels,
                           testCase.targetEdgeLabels, variant, domainFilterThreshold, "RI");
}

/// The refinement removes target nodes 2 and 3, half of the domain of pattern node 1, since they
/// have no neighbor with label 1. It removes nothing from the domain of pattern node 0.
Case halfRefinedCase() {
    return {.name = "half-refined",
            .pattern = IsomorphismTest::graphOf(2, {{0, 1}}),
            .target = IsomorphismTest::graphOf(6, {{0, 4}, {1, 5}, {2, 3}}),
            .semantics = Semantics::MONOMORPHISM,
            .patternNodeLabels = {1, 0},
            .targetNodeLabels = {0, 0, 0, 0, 1, 1}};
}

/// Runs the search and returns the matches.
std::vector<Match> matchesOf(const RISearchSetup &setup,
                             std::span<const index> patternNodeLabels = {},
                             std::span<const index> targetNodeLabels = {},
                             Semantics semantics = Semantics::MONOMORPHISM) {
    Aux::SignalHandler handler;
    std::vector<Match> matches;
    RIImpl(setup, patternNodeLabels, targetNodeLabels, semantics, handler,
           [&matches](std::span<const node> match) -> bool {
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

    for (const int variant : {2, 256})
        EXPECT_THAT(
            [&]() -> void {
                RI algo(pattern, target, Semantics::INDUCED, 0, static_cast<RI::Variant>(variant));
            },
            ThrowsMessage<std::invalid_argument>(HasSubstr("variant")))
            << "variant " << variant;
    EXPECT_NO_THROW(RI algo(pattern, target, Semantics::INDUCED, 0, GetParam()));
}

/// The threshold is a fraction of a domain. NaN fails every comparison, so it gets a case of its
/// own.
TEST_P(RIGTest, testDomainFilterThresholdOutsideTheUnitIntervalIsRejected) {
    const Graph pattern = IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}});
    const Graph target = IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}, {2, 0}});

    for (const double threshold : {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN()})
        EXPECT_THAT(
            [&]() -> void {
                RI algo(pattern, target, Semantics::INDUCED, 0, GetParam(), threshold);
            },
            ThrowsMessage<std::invalid_argument>(HasSubstr("threshold")))
            << "threshold " << threshold;

    for (const double threshold : {0.0, 1.0})
        EXPECT_NO_THROW(RI algo(pattern, target, Semantics::INDUCED, 0, GetParam(), threshold))
            << "threshold " << threshold;
}

TEST_P(RIGTest, testAgreesWithTheReference) {
    const RI::Variant variant = GetParam();
    for (const double threshold : {0.0, defaultDomainFilterThreshold, 1.0}) {
        SCOPED_TRACE(testing::Message() << "threshold " << threshold);
        const auto make =
            [variant, threshold](const Graph &pattern, const Graph &target, Semantics semantics,
                                 count maxMatches) -> std::unique_ptr<SubgraphIsomorphism> {
            return std::make_unique<RI>(pattern, target, semantics, maxMatches, variant, threshold);
        };

        IsomorphismTest::expectMatchesReference(make);
        IsomorphismTest::expectRespectsMatchCap(make);
        IsomorphismTest::expectCallbackStopsTheSearch(make);
    }
}

TEST_P(RIGTest, testOrderingInvariants) {
    for (const Case &testCase : IsomorphismTest::standardCases()) {
        const RISearchSetup setup = setupOf(testCase, GetParam());
        const SearchGraph &pattern = setup.patternGraph;
        const RIOrdering &ordering = setup.ordering;

        // An impossible instance gets no order.
        if (!setup.canMatch) {
            EXPECT_THAT(ordering.order, IsEmpty()) << "case: " << testCase.name;
            continue;
        }

        // An assertion, since the loop below indexes the pattern with these nodes.
        const auto nodes = testCase.pattern.nodeRange();
        ASSERT_THAT(ordering.order, UnorderedElementsAreArray(nodes.begin(), nodes.end()))
            << "case: " << testCase.name;
        ASSERT_THAT(ordering.parent, SizeIs(ordering.order.size())) << "case: " << testCase.name;

        for (index i = 0; i < ordering.order.size(); ++i) {
            const node pu = ordering.order[i];

            std::vector<index> earlierNeighbors;
            bool mutualWithFirst = false;
            for (index j = 0; j < i; ++j) {
                const node pj = ordering.order[j];
                const bool forward = pattern.hasEdge(pj, pu);
                const bool backward = pattern.hasEdge(pu, pj);
                if (!forward && !backward)
                    continue;
                if (earlierNeighbors.empty())
                    mutualWithFirst = pattern.isDirected() && forward && backward;
                earlierNeighbors.push_back(j);
            }

            // A `none` parent means exactly that no earlier position is adjacent.
            EXPECT_EQ(ordering.parent[i], earlierNeighbors.empty() ? none : earlierNeighbors[0])
                << "case: " << testCase.name << " at position " << i;

            // The candidates come from the slice along the parent arc, so only a directed mutual
            // pair keeps the parent among the checked neighbors.
            if (!earlierNeighbors.empty() && !mutualWithFirst)
                earlierNeighbors.erase(earlierNeighbors.begin());
            EXPECT_THAT(ordering.checkedNeighborsOf(i), ElementsAreArray(earlierNeighbors))
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
    const auto orderingOf = [variant](const Graph &pattern) -> RIOrdering {
        return setupOf(pattern, pattern, {}, {}, variant).ordering;
    };

    // The middle node of a path has the maximum degree and goes first.
    const RIOrdering path = orderingOf(IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}}));
    EXPECT_THAT(path.order, ElementsAre(1u, 0u, 2u));
    EXPECT_THAT(path.parent, ElementsAre(none, 0u, 0u));

    // Position 2 starts a new component.
    const RIOrdering split = orderingOf(IsomorphismTest::graphOf(4, {{0, 1}, {2, 3}}));
    EXPECT_THAT(split.order, ElementsAre(0u, 1u, 2u, 3u));
    EXPECT_THAT(split.parent, ElementsAre(none, 0u, none, 2u));

    // In a triangle, the last node checks only its arc to the second node. Its candidates come from
    // the slice of its parent, the first node.
    const RIOrdering triangle = orderingOf(IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}, {2, 0}}));
    EXPECT_THAT(triangle.order, ElementsAre(0u, 1u, 2u));
    EXPECT_THAT(triangle.parent, ElementsAre(none, 0u, 0u));
    EXPECT_THAT(triangle.checkedNeighborsOf(2), ElementsAre(1u));

    // The candidates come from one arc of a mutual pair, so the other arc to the parent is checked.
    const RIOrdering mutual =
        orderingOf(IsomorphismTest::graphOf(2, {{0, 1}, {1, 0}}, /* directed = */ true));
    EXPECT_THAT(mutual.parent, ElementsAre(none, 0u));
    EXPECT_THAT(mutual.checkedNeighborsOf(1), ElementsAre(0u));

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
        EXPECT_THAT(tied.ordering.order, ElementsAre(0u, 1u))
            << "with no domains a full tie goes to the smallest node id";
        return;
    }

    ASSERT_THAT(tied.domains.ofPatternNode[0], SizeIs(3));
    ASSERT_THAT(tied.domains.ofPatternNode[1], SizeIs(2));
    EXPECT_THAT(tied.ordering.order, ElementsAre(1u, 0u))
        << "the tie on all three counts must go to the more constrained node";
}

/// Unlike on the corpus, the domains prune on karate. ChibaNishizekiTriangleEdgeScore provides an
/// independent expected count.
TEST_P(RIGTest, testVariantsAgreeOnKarate) {
    Graph karate = METISGraphReader{}.read("input/karate.graph");
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

    ASSERT_THAT(setupOf(pattern, karate, {}, {}, RI::Variant::RI_DS).domains.ofPatternNode[0],
                SizeIs(Lt(karate.numberOfNodes())));

    RI plain(pattern, karate, Semantics::MONOMORPHISM, 0, RI::Variant::RI);
    plain.run();

    RI withDomains(pattern, karate, Semantics::MONOMORPHISM, 0, RI::Variant::RI_DS);
    withDomains.run();

    RI filteringEverywhere(pattern, karate, Semantics::MONOMORPHISM, 0, RI::Variant::RI_DS, 0.0);
    filteringEverywhere.run();

    EXPECT_THAT(withDomains.getMatches(), UnorderedElementsAreArray(plain.getMatches()))
        << "RI-DS is pure pruning and must not change the match set";
    EXPECT_THAT(filteringEverywhere.getMatches(), UnorderedElementsAreArray(plain.getMatches()))
        << "filtering the candidates through every domain must not change the match set";
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
    ASSERT_THAT(setup.ordering.parent, ElementsAre(none, 0u));
    ASSERT_TRUE(setup.domains.filtersCandidates[1])
        << "RI-DS no longer intersects the candidates of pattern node 1, so this case tests "
           "nothing";

    const std::vector<Match> expected = IsomorphismTest::referenceMatches(
        pattern, target, Semantics::MONOMORPHISM, patternNodeLabels, targetNodeLabels);
    ASSERT_THAT(expected, Not(IsEmpty())) << "a case with no matches would not exercise the search";

    RI algo(pattern, target, Semantics::MONOMORPHISM, 0, GetParam());
    algo.setNodeLabels(patternNodeLabels, targetNodeLabels);
    algo.run();

    EXPECT_THAT(algo.getMatches(), UnorderedElementsAreArray(expected));
}

/// The domains are set by hand, since a computed domain keeps every neighbor that could host the
/// pattern node. Pattern node 1 has node 0 as its parent, so its candidates come from the slice of
/// target node 5. That slice has a node below and a node above the domain.
TEST_F(RIDSGTest, testCandidatesAreIntersectedWithTheDomain) {
    const SearchGraph pattern(IsomorphismTest::graphOf(2, {{0, 1}}), /* buildMatrix = */ true);
    const RISearchSetup setup{
        .patternGraph = pattern,
        .targetGraph = SearchGraph(IsomorphismTest::graphOf(10, {{5, 1}, {5, 3}, {5, 7}, {5, 9}}),
                                   /* buildMatrix = */ false),
        .domains = {.ofPatternNode = {{5}, {3, 7}}, .filtersCandidates = {false, true}},
        .ordering = orderingFrom(pattern, {0, 1})};

    EXPECT_THAT(matchesOf(setup), UnorderedElementsAre(Match{5, 3}, Match{5, 7}));
}

/// The directed counterpart. Pattern node 1 points to its parent, so its candidates come from the
/// in-slice of target node 5. Domain node 4 is only an out-neighbor of node 5 and must not appear.
TEST_F(RIDSGTest, testInSliceCandidatesAreIntersectedWithTheDomain) {
    const SearchGraph pattern(IsomorphismTest::graphOf(2, {{1, 0}}, /* directed = */ true),
                              /* buildMatrix = */ true);
    const RISearchSetup setup{
        .patternGraph = pattern,
        .targetGraph =
            SearchGraph(IsomorphismTest::graphOf(10, {{1, 5}, {3, 5}, {7, 5}, {9, 5}, {5, 4}},
                                                 /* directed = */ true),
                        /* buildMatrix = */ false),
        .domains = {.ofPatternNode = {{5}, {3, 4, 7}}, .filtersCandidates = {false, true}},
        .ordering = orderingFrom(pattern, {0, 1})};

    EXPECT_THAT(matchesOf(setup), UnorderedElementsAre(Match{5, 3}, Match{5, 7}));
}

/// No domain of this case ends with a single node, so only the threshold decides which domains
/// filter the candidates.
TEST_F(RIDSGTest, testDomainFilterThresholdSelectsTheDomains) {
    const Case testCase = halfRefinedCase();
    const auto filtering = [&testCase](double threshold) -> std::vector<bool> {
        return setupOf(testCase, RI::Variant::RI_DS, threshold).domains.filtersCandidates;
    };

    ASSERT_THAT(setupOf(testCase, RI::Variant::RI_DS).domains.ofPatternNode,
                ElementsAre(ElementsAre(4u, 5u), ElementsAre(0u, 1u)));

    EXPECT_THAT(filtering(0.0), ElementsAre(true, true));
    EXPECT_THAT(filtering(0.5), ElementsAre(false, true));
    EXPECT_THAT(filtering(0.6), ElementsAre(false, false));
    EXPECT_THAT(filtering(1.0), ElementsAre(false, false));
}

/**
 * A singleton domain opens the matching order, and forward checking removes its target node from
 * all other domains. Pattern node 2 is isolated, so the refinement cannot relate it to nodes 0 and
 * 1. Only forward checking can then remove target node 3 from their domains.
 */
TEST_P(RIGTest, testSingletonDomainsAreRemovedAndComeFirst) {
    const auto expectSingletonsFirst = [](const RISearchSetup &prep, const char *name) -> void {
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
        EXPECT_THAT(prep.domains.ofPatternNode, IsEmpty())
            << "plain RI must build no domains at all";
        EXPECT_EQ(prep.ordering.order.front(), 0u)
            << "plain RI has no domains, so it starts at a maximum-degree node as it always did";
        return;
    }

    const std::vector<std::vector<node>> &domains = prep.domains.ofPatternNode;

    EXPECT_TRUE(prep.canMatch);
    EXPECT_THAT(domains[2], ElementsAre(unique));
    EXPECT_THAT(domains[0], ElementsAre(0u, 1u, 2u, 4u, 5u));
    EXPECT_THAT(domains[1], ElementsAre(0u, 1u, 2u, 4u, 5u))
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
        target.forNodes([&](node v) -> void {
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

    const RISearchSetup setup =
        setupOf(pattern, target, patternNodeLabels, targetNodeLabels, RI::Variant::RI_DS);
    const RIDomains &domains = setup.domains;

    EXPECT_TRUE(setup.canMatch);
    for (node center = 0; center < 4; ++center)
        EXPECT_THAT(domains.ofPatternNode[center], ElementsAre(center))
            << "center " << center << " should have been narrowed to one target node by the "
            << "cascade";
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

    EXPECT_FALSE(
        setupOf(pattern, target, patternNodeLabels, targetNodeLabels, RI::Variant::RI_DS).canMatch)
        << "one target node cannot host two pattern nodes at once";

    RI algo(pattern, target, Semantics::MONOMORPHISM, 0, GetParam());
    algo.setNodeLabels(patternNodeLabels, targetNodeLabels);
    algo.run();
    EXPECT_EQ(algo.numberOfMatches(), 0u);
}

/// Too few target nodes or too small target degrees rule out every match before the domains and the
/// order are computed.
TEST_P(RIGTest, testSizesAndDegreesRuleOutEveryMatch) {
    struct Impossible {
        const char *name;
        Graph pattern;
        Graph target;
    };

    const std::vector<Impossible> cases{
        {"more pattern nodes", IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}}),
         IsomorphismTest::graphOf(2, {{0, 1}})},
        {"larger degree", IsomorphismTest::graphOf(4, {{0, 1}, {0, 2}, {0, 3}}),
         IsomorphismTest::graphOf(4, {{0, 1}, {1, 2}, {2, 3}})},
        {"larger in-degree",
         IsomorphismTest::graphOf(4, {{1, 0}, {2, 0}, {3, 0}}, /* directed = */ true),
         IsomorphismTest::graphOf(4, {{0, 1}, {0, 2}, {0, 3}}, /* directed = */ true)}};

    for (const Impossible &impossible : cases) {
        const RISearchSetup setup =
            setupOf(impossible.pattern, impossible.target, {}, {}, GetParam());
        EXPECT_FALSE(setup.canMatch) << impossible.name;
        EXPECT_THAT(setup.domains.ofPatternNode, IsEmpty()) << impossible.name;
        EXPECT_THAT(setup.ordering.order, IsEmpty()) << impossible.name;
    }

    // The empty mapping matches an empty pattern.
    EXPECT_TRUE(setupOf(Graph(0), Graph(0), {}, {}, GetParam()).canMatch);
}

/// Drives RIImpl one level at a time, as the ParallelRI workers do. A difference to run() means
/// that a state lacks something that run() carries from one depth to the next.
TEST_P(RIGTest, testExpandAgreesWithRun) {
    for (const Case &testCase : IsomorphismTest::standardCases()) {
        const RISearchSetup setup = setupOf(testCase, GetParam());
        // ParallelRI checks this before it calls expand().
        if (!setup.canMatch)
            continue;

        std::vector<Match> viaRun = matchesOf(setup, testCase.patternNodeLabels,
                                              testCase.targetNodeLabels, testCase.semantics);

        Aux::SignalHandler handler;
        std::vector<Match> viaExpand;
        RIImpl expander(setup, testCase.patternNodeLabels, testCase.targetNodeLabels,
                        testCase.semantics, handler,
                        [&viaExpand](std::span<const node> match) -> bool {
                            viaExpand.emplace_back(match.begin(), match.end());
                            return true;
                        });

        std::vector<RIState> pending{expander.rootState()};
        while (!pending.empty()) {
            // Moved out first, since expand() appends the children to `pending`.
            const RIState state = std::move(pending.back());
            pending.pop_back();
            EXPECT_TRUE(expander.expand(state, pending)) << "case: " << testCase.name;
        }

        // Sorted, as in expectMatchesReference(): UnorderedElementsAreArray would compare every
        // pair of matches, and some random cases have over a thousand matches.
        std::ranges::sort(viaRun);
        EXPECT_THAT(viaExpand, WhenSorted(ElementsAreArray(viaRun))) << "case: " << testCase.name;
    }
}

/// The search checks for an interrupt before every report, so no match follows the interrupt.
TEST_P(RIGTest, testInterruptStopsTheSearch) {
    const Graph pattern = IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}});
    const Graph cycle =
        IsomorphismTest::graphOf(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}});

    RI algo(pattern, cycle, Semantics::MONOMORPHISM, 0, GetParam());
    count delivered = 0;
    algo.setCallback([&](std::span<const node>) -> void {
        if (++delivered == 5)
            GlobalState::setReceivedSIGINT(true);
    });

    EXPECT_THROW(algo.run(), Aux::SignalHandler::InterruptException);
    GlobalState::setReceivedSIGINT(false);

    EXPECT_EQ(delivered, 5u);
    EXPECT_FALSE(algo.hasFinished());
}

/// One combination of the random cases: the variant, directed graphs, the semantics, node labels
/// and the domain filter threshold. Under RI-DS, threshold 0 filters through every domain and 1
/// only through the single-node ones.
using RandomCombination = std::tuple<RI::Variant, bool, Semantics, bool, double>;

class RIRandomGTest : public testing::TestWithParam<RandomCombination> {};

/// A test name such as `RIDS_directed_induced_labelled_threshold80`, which --gtest_filter can
/// select.
std::string combinationName(const testing::TestParamInfo<RandomCombination> &info) {
    const auto [variant, directed, semantics, nodeLabelled, threshold] = info.param;
    return std::string(variant == RI::Variant::RI ? "RI" : "RIDS")
           + (directed ? "_directed" : "_undirected")
           + (semantics == Semantics::INDUCED ? "_induced" : "_monomorphism")
           + (nodeLabelled ? "_labelled" : "_unlabelled") + "_threshold"
           + std::to_string(std::lround(threshold * 100));
}

INSTANTIATE_TEST_SUITE_P(
    Combinations, RIRandomGTest,
    testing::Combine(testing::Values(RI::Variant::RI, RI::Variant::RI_DS), testing::Bool(),
                     testing::Values(Semantics::MONOMORPHISM, Semantics::INDUCED), testing::Bool(),
                     testing::Values(0.0, defaultDomainFilterThreshold, 1.0)),
    combinationName);

TEST_P(RIRandomGTest, testMatchesReference) {
    const auto [variant, directed, semantics, nodeLabelled, threshold] = GetParam();

    // The reference enumerates all 6^4 mappings per trial, which keeps the sizes small.
    constexpr count patternNodes = 4;
    constexpr count targetNodes = 6;
    constexpr int trials = 8;

    Aux::Random::setSeed(1701, false);

    for (int trial = 0; trial < trials; ++trial) {
        const Graph pattern = ErdosRenyiGenerator(patternNodes, 0.5, directed).generate();
        const Graph target = ErdosRenyiGenerator(targetNodes, 0.4, directed).generate();

        std::vector<index> patternNodeLabels;
        std::vector<index> targetNodeLabels;
        if (nodeLabelled) {
            patternNodeLabels = randomNodeLabels(pattern.upperNodeIdBound());
            targetNodeLabels = randomNodeLabels(target.upperNodeIdBound());
        }

        const std::vector<Match> expected = IsomorphismTest::referenceMatches(
            pattern, target, semantics, patternNodeLabels, targetNodeLabels);

        RI algo(pattern, target, semantics, 0, variant, threshold);
        if (nodeLabelled)
            algo.setNodeLabels(patternNodeLabels, targetNodeLabels);
        algo.run();

        EXPECT_THAT(algo.getMatches(), UnorderedElementsAreArray(expected)) << "trial " << trial;
    }
}

} // namespace
} // namespace NetworKit
