/*  BFSTraversalGTest.cpp
 *
 *  Created on: 25.10.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/BFS.hpp>

#include "TraversalTestHelper.hpp"

namespace NetworKit {
namespace {

using ::testing::ElementsAre;
using ::testing::FieldsAre;
using ::testing::IsEmpty;

template <class TestT>
using GenericBFSTraversalGTest = GenericTraversalGTest<TestT>;

TYPED_TEST_SUITE(GenericBFSTraversalGTest, TraversalTestTypes,
                 /*Comma needed for variadic macro.*/);

TYPED_TEST(GenericBFSTraversalGTest, testBFSfromTypedGraphs) {
    using NodeT = typename TestFixture::NodeT;

    const auto G = this->unweightedGraph();

    std::vector<NodeT> sourceSequence;
    std::vector<count> sourceDistances;
    Traversal::BFSfrom(G, NodeT{0}, [&](NodeT u, count dist) {
        sourceSequence.push_back(u);
        sourceDistances.push_back(dist);
    });

    EXPECT_THAT(sourceSequence,
                ElementsAre(NodeT{0}, NodeT{1}, NodeT{2}, NodeT{3}, NodeT{4}, NodeT{5}, NodeT{6}));
    EXPECT_THAT(sourceDistances, ElementsAre(0, 1, 1, 2, 2, 3, 3));

    const std::vector<NodeT> sources{NodeT{0}, NodeT{6}};
    std::vector<NodeT> rangeSequence;
    Traversal::BFSfrom(G, sources.begin(), sources.end(),
                       [&](NodeT u) { rangeSequence.push_back(u); });

    EXPECT_THAT(rangeSequence,
                ElementsAre(NodeT{0}, NodeT{6}, NodeT{1}, NodeT{2}, NodeT{4}, NodeT{3}, NodeT{5}));
}

TYPED_TEST(GenericBFSTraversalGTest, testBFSfromEmptySourceRange) {
    using NodeT = typename TestFixture::NodeT;

    const auto G = this->unweightedGraph();
    const std::vector<NodeT> sources;

    count calls = 0;
    Traversal::BFSfrom(G, sources.begin(), sources.end(), [&](NodeT) { ++calls; });

    EXPECT_EQ(calls, 0);
}

TYPED_TEST(GenericBFSTraversalGTest, testBFSfromDuplicateSources) {
    using NodeT = typename TestFixture::NodeT;

    const auto G = this->pathGraph();
    const std::vector<NodeT> sources{NodeT{0}, NodeT{0}, NodeT{2}};

    std::vector<NodeT> sequence;
    std::vector<count> distances;
    Traversal::BFSfrom(G, sources.begin(), sources.end(), [&](NodeT u, count dist) {
        sequence.push_back(u);
        distances.push_back(dist);
    });

    EXPECT_THAT(sequence, ElementsAre(NodeT{0}, NodeT{2}, NodeT{1}, NodeT{3}, NodeT{4}));
    EXPECT_THAT(distances, ElementsAre(0, 0, 1, 1, 2));
}

TYPED_TEST(GenericBFSTraversalGTest, testBFSfromInvalidSources) {
    using NodeT = typename TestFixture::NodeT;

    auto G = this->pathGraph();
    const std::vector<NodeT> outOfRangeSources{NodeT{0}, static_cast<NodeT>(G.upperNodeIdBound())};
    EXPECT_THROW(
        Traversal::BFSfrom(G, outOfRangeSources.begin(), outOfRangeSources.end(), [](NodeT) {}),
        std::runtime_error);

    G.removeNode(NodeT{2});
    const std::vector<NodeT> removedNodeSources{NodeT{0}, NodeT{2}};
    EXPECT_THROW(
        Traversal::BFSfrom(G, removedNodeSources.begin(), removedNodeSources.end(), [](NodeT) {}),
        std::runtime_error);

    if constexpr (std::is_signed_v<NodeT>) {
        const std::vector<NodeT> negativeSources{NodeT{0}, NodeT{-1}};
        EXPECT_THROW(
            Traversal::BFSfrom(G, negativeSources.begin(), negativeSources.end(), [](NodeT) {}),
            std::runtime_error);
    }
}

TYPED_TEST(GenericBFSTraversalGTest, testBFSfromMultiSourceDistances) {
    using NodeT = typename TestFixture::NodeT;

    const auto G = this->pathGraph();
    const std::vector<NodeT> sources{NodeT{0}, NodeT{4}};

    std::vector<NodeT> sequence;
    std::vector<count> distances;
    Traversal::BFSfrom(G, sources.begin(), sources.end(), [&](NodeT u, count dist) {
        sequence.push_back(u);
        distances.push_back(dist);
    });

    EXPECT_THAT(sequence, ElementsAre(NodeT{0}, NodeT{4}, NodeT{1}, NodeT{3}, NodeT{2}));
    EXPECT_THAT(distances, ElementsAre(0, 0, 1, 1, 2));
}

TYPED_TEST(GenericBFSTraversalGTest, testBFSfromDirectedGraph) {
    using NodeT = typename TestFixture::NodeT;

    const auto G = this->directedGraph();
    std::vector<NodeT> sequence;
    Traversal::BFSfrom(G, NodeT{0}, [&](NodeT u) { sequence.push_back(u); });

    EXPECT_THAT(sequence, ElementsAre(NodeT{0}, NodeT{1}, NodeT{2}, NodeT{4}));
}

TYPED_TEST(GenericBFSTraversalGTest, testBFSEdgesFromTypedGraphs) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    auto G = this->weightedGraph();
    G.indexEdges();
    std::vector<std::tuple<NodeT, NodeT, EdgeWeightT, edgeid>> edgeSequence;

    Traversal::BFSEdgesFrom(G, NodeT{0}, [&](NodeT u, NodeT v, EdgeWeightT w, edgeid eid) {
        edgeSequence.emplace_back(u, v, w, eid);
    });

    EXPECT_THAT(
        edgeSequence,
        ElementsAre(FieldsAre(NodeT{0}, NodeT{1}, EdgeWeightT{10}, G.edgeId(NodeT{0}, NodeT{1})),
                    FieldsAre(NodeT{0}, NodeT{2}, EdgeWeightT{20}, G.edgeId(NodeT{0}, NodeT{2})),
                    FieldsAre(NodeT{1}, NodeT{3}, EdgeWeightT{30}, G.edgeId(NodeT{1}, NodeT{3})),
                    FieldsAre(NodeT{2}, NodeT{4}, EdgeWeightT{40}, G.edgeId(NodeT{2}, NodeT{4})),
                    FieldsAre(NodeT{3}, NodeT{5}, EdgeWeightT{50}, G.edgeId(NodeT{3}, NodeT{5})),
                    FieldsAre(NodeT{4}, NodeT{6}, EdgeWeightT{60}, G.edgeId(NodeT{4}, NodeT{6}))));
}

TYPED_TEST(GenericBFSTraversalGTest, testBFSEdgesFromCornerCases) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    GraphT G(5, true, true);
    G.addEdge(NodeT{0}, NodeT{0}, EdgeWeightT{10});
    G.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{20});
    G.addEdge(NodeT{0}, NodeT{2}, EdgeWeightT{30});
    G.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{40});
    G.addEdge(NodeT{3}, NodeT{0}, EdgeWeightT{50});
    G.indexEdges();

    std::vector<std::tuple<NodeT, NodeT, EdgeWeightT, edgeid>> edgeSequence;
    Traversal::BFSEdgesFrom(G, NodeT{0}, [&](NodeT u, NodeT v, EdgeWeightT w, edgeid eid) {
        edgeSequence.emplace_back(u, v, w, eid);
    });

    EXPECT_THAT(
        edgeSequence,
        ElementsAre(FieldsAre(NodeT{0}, NodeT{1}, EdgeWeightT{20}, G.edgeId(NodeT{0}, NodeT{1})),
                    FieldsAre(NodeT{0}, NodeT{2}, EdgeWeightT{30}, G.edgeId(NodeT{0}, NodeT{2}))));

    edgeSequence.clear();
    Traversal::BFSEdgesFrom(G, NodeT{4}, [&](NodeT u, NodeT v, EdgeWeightT w, edgeid eid) {
        edgeSequence.emplace_back(u, v, w, eid);
    });

    EXPECT_THAT(edgeSequence, IsEmpty());
}

TYPED_TEST(GenericBFSTraversalGTest, testBFSEdgesFromInvalidSource) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    const auto G = this->weightedGraph();
    EXPECT_THROW(Traversal::BFSEdgesFrom(G, static_cast<NodeT>(G.upperNodeIdBound()),
                                         [](NodeT, NodeT, EdgeWeightT, edgeid) {}),
                 std::runtime_error);
}

} // namespace
} // namespace NetworKit
