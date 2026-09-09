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

/// The two snapshots RI::run() builds, so a test can drive RIImpl exactly the way the driver does.
struct Snapshot {
    SearchGraph pattern;
    SearchGraph target;

    explicit Snapshot(const Case &testCase)
        : pattern(testCase.pattern, /* buildMatrix = */ true, testCase.patternEdgeLabels),
          target(testCase.target, /* buildMatrix = */ false, testCase.targetEdgeLabels) {}

    /// Whether RI::run() would refuse this case rather than search it - see
    /// SearchGraph::collapsedLabelledEdges().
    bool refused() const {
        return pattern.collapsedLabelledEdges() || target.collapsedLabelledEdges();
    }
};

/// The two preprocessing results, computed in the order both drivers compute them: the domains
/// first, then the order read off them. Doing it any other way would not test what ships.
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

/// Where @a pu sits in the matching order, or `none` if it is not in it.
index positionOf(const RIImpl::Ordering &ordering, node pu) {
    const auto found = std::find(ordering.order.begin(), ordering.order.end(), pu);
    return found == ordering.order.end() ? none
                                         : static_cast<index>(found - ordering.order.begin());
}

/// A label per node id, with `none` mixed in so both sides carry wildcards.
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

/**
 * Everything here is parameterised over the variant, so neither is tested less than the other.
 *
 * The two no longer produce the same matching order: RI-Ds puts singleton-domain nodes first and
 * settles a full tie by domain size, so its order depends on the target. Running a test twice is
 * therefore no longer evidence that the orders agree - where they must differ, the expectation
 * branches on the parameter, and where the domains are empty under plain RI the test says so
 * rather than asserting nothing. What running twice does still prove is the thing that matters
 * most: both variants must return exactly the same *matches*, because domains are pure pruning.
 */
class RIGTest : public testing::TestWithParam<RI::Variant> {};

INSTANTIATE_TEST_SUITE_P(Variants, RIGTest, testing::Values(RI::Variant::RI, RI::Variant::RI_DS));

// -------------------------------------------------------------------------------------------
// The three assertions the shared harness offers
// -------------------------------------------------------------------------------------------

TEST_P(RIGTest, testMatchesReference) {
    const RI::Variant variant = GetParam();
    const auto make = [variant](const Graph &pattern, const Graph &target, Semantics semantics,
                                count maxMatches) {
        return std::unique_ptr<SubgraphIsomorphism>(
            new RI(pattern, target, variant, semantics, maxMatches));
    };

    IsomorphismTest::expectMatchesReference(make);
}

TEST_P(RIGTest, testRespectsMatchCap) {
    const RI::Variant variant = GetParam();
    const auto make = [variant](const Graph &pattern, const Graph &target, Semantics semantics,
                                count maxMatches) {
        return std::unique_ptr<SubgraphIsomorphism>(
            new RI(pattern, target, variant, semantics, maxMatches));
    };

    IsomorphismTest::expectRespectsMatchCap(make);
}

TEST_P(RIGTest, testCallbackFormsAgree) {
    const RI::Variant variant = GetParam();
    const auto make = [variant](const Graph &pattern, const Graph &target, Semantics semantics,
                                count maxMatches) {
        return std::unique_ptr<SubgraphIsomorphism>(
            new RI(pattern, target, variant, semantics, maxMatches));
    };

    IsomorphismTest::expectCallbackFormsAgree(make);
}

// -------------------------------------------------------------------------------------------
// The ordering, which agreement with the reference judges only by luck on a 4-node corpus
// -------------------------------------------------------------------------------------------

/**
 * Structural invariants of the matching order, checked directly rather than through a search.
 *
 * This is what catches the whole "the ordering is subtly wrong" class: a removed id in the order,
 * a node listed twice, a parent pointing forwards, a parent that is not actually adjacent, or a
 * `none` parent at a position that does have an earlier neighbour.
 */
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

        if (!ordering.parent.empty())
            EXPECT_EQ(ordering.parent[0], none) << "case: " << testCase.name;

        for (index i = 0; i < ordering.order.size(); ++i) {
            const node pu = ordering.order[i];

            // A parent must be an *earlier* position that really is adjacent, in either direction.
            const index parentPos = ordering.parent[i];
            if (parentPos != none) {
                ASSERT_LT(parentPos, i) << "case: " << testCase.name;
                const node pp = ordering.order[parentPos];
                EXPECT_TRUE(snapshot.pattern.hasEdge(pp, pu) || snapshot.pattern.hasEdge(pu, pp))
                    << "case: " << testCase.name << " - parent is not adjacent";
            }

            // And `none` means exactly "this position starts a component", never "I gave up".
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
 * Two orders small enough to work out by hand, plus the pair that proves the score is three-level.
 *
 * The eight-node graph is built in the shape of the paper's worked example. At the step that
 * decides between nodes 5 and 0 the two tie on the first term (one arc into the order each) and on
 * the third (one untouched neighbour each), and differ only on the middle one: 5 reaches the
 * already-ordered node 1 through the still-unordered node 2, while 0 reaches nothing. So 5 coming
 * before 0 is the only direct evidence that the middle term is computed at all - a two-level score
 * would order them the other way round, by node id.
 */
TEST_P(RIGTest, testOrderingHandTraced) {
    const RI::Variant variant = GetParam();
    const auto orderingOf = [variant](const Graph &pattern) {
        const SearchGraph patternGraph(pattern, /* buildMatrix = */ true);
        // Under plain RI the target is never consulted; under RI-Ds it reaches the order through
        // the domain sizes, and passing the pattern as its own target keeps these three graphs
        // small enough to work out by hand.
        const SearchGraph targetGraph(pattern, /* buildMatrix = */ false);
        return Preprocessed(patternGraph, targetGraph, {}, {}, variant).ordering;
    };

    // A path 0-1-2. The middle node is the only one of degree two, so it goes first, and both ends
    // then hang off position 0.
    const RIImpl::Ordering path = orderingOf(IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}}));
    EXPECT_EQ(path.order, (std::vector<node>{1, 0, 2}));
    EXPECT_EQ(path.parent, (std::vector<index>{none, 0, 0}));

    // Two components. Position 2 restarts, which is what the `none` parent is for.
    const RIImpl::Ordering split = orderingOf(IsomorphismTest::graphOf(4, {{0, 1}, {2, 3}}));
    EXPECT_EQ(split.order, (std::vector<node>{0, 1, 2, 3}));
    EXPECT_EQ(split.parent, (std::vector<index>{none, 0, none, 2}));

    const RIImpl::Ordering worked = orderingOf(IsomorphismTest::graphOf(
        8, {{4, 1}, {4, 5}, {4, 0}, {4, 3}, {1, 2}, {1, 6}, {5, 2}, {5, 7}, {0, 7}, {3, 6}}));
    EXPECT_EQ(worked.order[0], 4u) << "the first node must be the unique maximum-degree node";
    EXPECT_LT(positionOf(worked, 5), positionOf(worked, 0))
        << "node 5 beats node 0 only on the two-hop term";
}

// -------------------------------------------------------------------------------------------
// RI-Ds, on a target big enough for a domain to actually remove something
// -------------------------------------------------------------------------------------------

/**
 * The corpus is far too small for a domain to prune, so RI-Ds is effectively untested by it.
 *
 * Domains are pure pruning, so any divergence between the variants is an unsound domain - the
 * refinement sweep removed a target node a real match needs. The absolute count is pinned against
 * ChibaNishizeki's per-edge triangle counts, which is an independent check from a part of
 * NetworKit that has nothing to do with this module.
 */
TEST_P(RIGTest, testVariantsAgreeOnKarate) {
    METISGraphReader reader;
    Graph karate = reader.read("input/karate.graph");
    karate.indexEdges();

    ChibaNishizekiTriangleEdgeScore triangleScore(karate);
    triangleScore.run();
    const std::vector<count> perEdge = triangleScore.scores();
    const count edgeSum = std::accumulate(perEdge.begin(), perEdge.end(), count{0});
    ASSERT_GT(edgeSum, 0u) << "karate should contain triangles";

    // Each triangle is counted on three edges and found once per way of labelling its three nodes,
    // so the number of matches is 6 * (edgeSum / 3) = 2 * edgeSum.
    const count expected = 2 * edgeSum;

    const Graph pattern = IsomorphismTest::graphOf(3, {{0, 1}, {1, 2}, {2, 0}});

    // The domains really do remove something here, which is the whole reason this case exists:
    // karate has nodes of degree one, and none of those can host a triangle node of degree two.
    // On the 4-node corpus above no domain ever loses an entry.
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
 * A target built so that the refinement sweep is selective enough for RI-Ds to intersect with.
 *
 * RI-Ds only applies a domain to a target slice when the sweep removed most of it, because below
 * that the domain re-rejects what the cheap rules reject anyway. On every graph in the corpus the
 * sweep removes nothing, so without this case the intersecting path would never run in a test at
 * all.
 *
 * The construction: ten nodes of class 0 and ten of class 1, but exactly one edge joins the two
 * classes. So a class-0 node can host the pattern's class-0 end only if it is that one node, and
 * the sweep drops nine of the ten. The pattern's class-0 node is given the higher id so that it
 * lands at the second position in the matching order, which is the parented one - the position
 * whose candidates come from a slice.
 *
 * The yield is asserted rather than assumed, so that an edit which quietly makes the graph less
 * selective fails here instead of silently stopping testing the path.
 */
TEST_P(RIGTest, testSelectiveDomainsDoNotChangeTheMatchSet) {
    constexpr count perClass = 10;

    Graph target(2 * perClass);
    // Class-0 nodes 1..9 are joined to each other, so they have degree but no class-1 neighbour.
    for (node u = 1; u + 1 < perClass; ++u)
        target.addEdge(u, u + 1);
    // The single edge across the classes.
    target.addEdge(0, perClass);

    std::vector<index> targetNodeLabels(2 * perClass);
    for (node v = 0; v < 2 * perClass; ++v)
        targetNodeLabels[v] = v < perClass ? 0 : 1;

    // Pattern node 1 carries class 0, so the order puts it second - at the parented position.
    const Graph pattern = IsomorphismTest::graphOf(2, {{0, 1}});
    const std::vector<index> patternNodeLabels{1, 0};

    // The precondition the whole case rests on: of the class-0 nodes that could host the pattern's
    // class-0 node on degree alone, at most a fifth survive the sweep's "must reach the other
    // position's domain" test. Anything above a fifth and RI-Ds stops intersecting.
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

// -------------------------------------------------------------------------------------------
// Forward checking and the domain-size ordering, which are RI-DS-SI-FC's whole content
// -------------------------------------------------------------------------------------------

/**
 * A singleton domain names the only target node its pattern node can take, so injectivity forbids
 * that node to every other pattern node. This asserts the removal happens.
 *
 * The construction keeps the refinement sweep out of it. Pattern node 2 is **isolated**, so it is
 * adjacent to neither of the others and the sweep - which only ever relates adjacent pattern nodes
 * - has no way to reach them from it. The unique label makes its domain a singleton, and forward
 * checking is then the only thing in the pipeline that can take that target node away from nodes 0
 * and 1. The assertion below that target node 2 survives in a domain is what pins that down: it is
 * a neighbour of target node 3, which is exactly the reason the sweep had to keep 3.
 */
TEST_P(RIGTest, testForwardCheckingRemovesSingletonImages) {
    constexpr node unique = 3;
    constexpr index rareLabel = 7;

    const Graph target = IsomorphismTest::graphOf(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});
    std::vector<index> targetNodeLabels(6, 0);
    targetNodeLabels[unique] = rareLabel;

    // Node 2 is isolated and carries the rare label; nodes 0 and 1 are wildcards, so nothing but
    // forward checking can keep them off target node 3.
    const Graph pattern = IsomorphismTest::graphOf(3, {{0, 1}});
    const std::vector<index> patternNodeLabels{none, none, rareLabel};

    const SearchGraph patternGraph(pattern, /* buildMatrix = */ true);
    const SearchGraph targetGraph(target, /* buildMatrix = */ false);
    const RIImpl::Domains domains = RIImpl::computeDomains(
        patternGraph, targetGraph, patternNodeLabels, targetNodeLabels, GetParam());

    if (GetParam() == RI::Variant::RI) {
        EXPECT_TRUE(domains.ofPatternNode.empty()) << "plain RI must build no domains at all";
        return;
    }

    EXPECT_FALSE(domains.anyEmpty);
    EXPECT_EQ(domains.ofPatternNode[2], (std::vector<node>{unique}));
    EXPECT_EQ(domains.ofPatternNode[0], (std::vector<node>{0, 1, 2, 4, 5}));
    EXPECT_EQ(domains.ofPatternNode[1], (std::vector<node>{0, 1, 2, 4, 5}));
    EXPECT_TRUE(std::find(domains.ofPatternNode[1].begin(), domains.ofPatternNode[1].end(), 2u)
                != domains.ofPatternNode[1].end())
        << "target node 2 is a neighbour of target node 3, which is why the sweep kept 3 - without "
           "this the removal above could be the sweep's rather than forward checking's";
}

/**
 * Removing a singleton's target node can leave some other domain with one entry, and that node's
 * own image is then decided too. One call to computeDomains() has to drain the whole cascade.
 *
 * Four pattern stars whose centres have degrees 4, 3, 2 and 1 sit over a target with exactly one
 * label-0 node of degree 4, two of degree 3 or more, three of degree 2 or more and four of degree 1
 * or more. Degree domination alone therefore gives the four centres strictly nested domains of
 * sizes 1, 2, 3 and 4, and the sweep cannot touch them because every one of those target nodes has
 * a label-9 neighbour that the leaves' domains still hold. So the three larger centres are
 * singletons at the end only if the removals chained: 1 forces 2, 2 forces 3, 3 forces 4.
 */
TEST_P(RIGTest, testForwardCheckingReachesAFixpoint) {
    constexpr index leafLabel = 9;

    // Target nodes 0..3 carry label 0 and have degrees 4, 3, 2 and 1; 4..13 are label-9 pendants.
    const Graph target = IsomorphismTest::graphOf(
        14, {{0, 4}, {0, 11}, {0, 12}, {0, 13}, {1, 5}, {1, 9}, {1, 10}, {2, 6}, {2, 8}, {3, 7}});
    std::vector<index> targetNodeLabels(14, leafLabel);
    for (node v = 0; v < 4; ++v)
        targetNodeLabels[v] = 0;

    // The nesting the cascade rests on, re-derived here from the target rather than assumed.
    for (const count degree : {4u, 3u, 2u, 1u}) {
        count labelZeroNodesOfThatDegree = 0;
        target.forNodes([&](node v) {
            if (targetNodeLabels[v] == 0 && target.degree(v) >= degree)
                ++labelZeroNodesOfThatDegree;
        });
        ASSERT_EQ(labelZeroNodesOfThatDegree, 5 - degree);
    }

    // Centres 0..3 with 4, 3, 2 and 1 leaves. The centres come first so that the sweep, which is a
    // single pass in pattern node id order, sees the leaves' domains before anything refined them.
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

/**
 * Two pattern nodes that can only go on the same target node make the instance impossible, and
 * forward checking is what notices before the search starts.
 *
 * Plain RI finds zero matches here by searching, without any domain machinery at all, which is what
 * makes this evidence of a correct early exit rather than of over-pruning: if RI-Ds were throwing
 * real matches away, plain RI would have found them.
 */
TEST_P(RIGTest, testForwardCheckingRejectsImpossibleInstances) {
    constexpr index rareLabel = 7;

    const Graph target = IsomorphismTest::graphOf(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});
    std::vector<index> targetNodeLabels(6, 0);
    targetNodeLabels[2] = rareLabel;

    // Two components, so the sweep never relates nodes 0 and 2 - both of which need target node 2.
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

/**
 * Section 4.1: a pattern node whose domain holds one target node is mapped before any node whose
 * domain is larger, because its image is already decided and mapping it first constrains
 * everything that follows for free.
 *
 * Checked over the whole corpus and over a labelled case built to have a singleton, because the
 * rule is a hard filter on the order and so has to hold at every position, not just the first.
 */
TEST_P(RIGTest, testSingletonDomainsComeFirst) {
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

    // The corpus is unlabelled, so a singleton only turns up when one is engineered. This is the
    // construction from testForwardCheckingRemovesSingletonImages: pattern node 2 is the singleton
    // and must be ordered first, ahead of the connected component that would otherwise start.
    const Graph target = IsomorphismTest::graphOf(6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});
    std::vector<index> targetNodeLabels(6, 0);
    targetNodeLabels[3] = 7;
    const Graph pattern = IsomorphismTest::graphOf(3, {{0, 1}});
    const std::vector<index> patternNodeLabels{none, none, 7};

    const SearchGraph patternGraph(pattern, /* buildMatrix = */ true);
    const SearchGraph targetGraph(target, /* buildMatrix = */ false);
    expectSingletonsFirst(patternGraph, targetGraph, patternNodeLabels, targetNodeLabels,
                          GetParam(), "engineered-singleton");

    const Preprocessed prep(patternGraph, targetGraph, patternNodeLabels, targetNodeLabels,
                            GetParam());
    if (GetParam() == RI::Variant::RI_DS)
        EXPECT_EQ(prep.ordering.order.front(), 2u)
            << "the singleton-domain node must open the order under RI-Ds";
    else
        EXPECT_EQ(prep.ordering.order.front(), 0u)
            << "plain RI has no domains, so it starts at a maximum-degree node as it always did";
}

/**
 * Section 4.2.1: a candidate that ties on the whole triple is settled by the smaller domain.
 *
 * Two isolated pattern nodes tie on all three terms - no arcs into the order, no two-hop reach, no
 * untouched neighbours - so under plain RI the order is decided by node id alone and node 0 wins.
 * Giving node 1 the rarer label leaves it the smaller domain without making it a singleton, so
 * flipping the order is the only thing the domain-size key can be doing. This is the analogue of
 * the "node 5 beats node 0 only on the two-hop term" assertion above, one level further down.
 */
TEST_P(RIGTest, testDomainSizeBreaksTiesInTheOrdering) {
    // No edges anywhere, so every node's triple is (0, 0, 0) and only the labels differ.
    const Graph target(5);
    const std::vector<index> targetNodeLabels{1, 1, 1, 2, 2};

    const Graph pattern(2);
    const std::vector<index> patternNodeLabels{1, 2};

    const SearchGraph patternGraph(pattern, /* buildMatrix = */ true);
    const SearchGraph targetGraph(target, /* buildMatrix = */ false);
    const Preprocessed prep(patternGraph, targetGraph, patternNodeLabels, targetNodeLabels,
                            GetParam());

    if (GetParam() == RI::Variant::RI) {
        EXPECT_EQ(prep.ordering.order, (std::vector<node>{0, 1}))
            << "with no domains a full tie goes to the smallest node id - erratum 6";
        return;
    }

    // Three target nodes can host pattern node 0, two can host pattern node 1. Neither is a
    // singleton, so the singletons-first rule is not what is being measured here.
    ASSERT_EQ(prep.domains.ofPatternNode[0].size(), 3u);
    ASSERT_EQ(prep.domains.ofPatternNode[1].size(), 2u);
    EXPECT_EQ(prep.ordering.order, (std::vector<node>{1, 0}))
        << "the tie on all three counts must go to the more constrained node";
}

// -------------------------------------------------------------------------------------------
// expand(), which nothing else exercises while ParallelRIImpl is stubbed
// -------------------------------------------------------------------------------------------

/**
 * Drive RIImpl a level at a time, in the shape ParallelRIImpl::workerLoop() will use, and assert
 * the match set is the recursion's.
 *
 * Nothing is carried on the C++ call stack here: every state on the pending list has to be
 * self-contained, which is exactly the property a stolen state needs. A disagreement means
 * something the recursion keeps on its own stack was never written into the State.
 */
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

        std::vector<RIImpl::State> pending(1);
        pending.front().mapping.assign(ordering.order.size(), none);
        while (!pending.empty()) {
            RIImpl::State state = std::move(pending.back());
            pending.pop_back();

            std::vector<RIImpl::State> children;
            // No cap is set, so the only reason expand() could ask to stop is a bug.
            EXPECT_TRUE(expander.expand(state, children)) << "case: " << testCase.name;

            for (RIImpl::State &child : children)
                pending.push_back(std::move(child));
        }

        IsomorphismTest::sortMatches(viaRun);
        IsomorphismTest::sortMatches(viaExpand);
        EXPECT_EQ(viaExpand, viaRun) << "case: " << testCase.name;
    }
}

// -------------------------------------------------------------------------------------------
// Randomised cross-check, for the ordering or parent bug the hand-built corpus happens to miss
// -------------------------------------------------------------------------------------------

TEST_P(RIGTest, testMatchesReferenceOnRandomGraphs) {
    const RI::Variant variant = GetParam();

    // The reference enumerates every injective mapping, so the sizes here are what keeps this in
    // the millisecond range: 6 target nodes to the power of 4 pattern nodes.
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
