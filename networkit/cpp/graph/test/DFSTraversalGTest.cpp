/*  DFSTraversalGTest.cpp
 *
 *  Created on: 25.10.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/DFS.hpp>

#include "TraversalTestHelper.hpp"

namespace NetworKit {
namespace {

using ::testing::ElementsAre;
using ::testing::FieldsAre;
using ::testing::IsEmpty;
using ::testing::Pair;

template <class TestT>
using GenericDFSTraversalGTest = GenericTraversalGTest<TestT>;

TYPED_TEST_SUITE(GenericDFSTraversalGTest, TraversalTestTypes,
                 /*Comma needed for variadic macro.*/);

TYPED_TEST(GenericDFSTraversalGTest, testDFSfromTypedGraphs) {
    using NodeT = typename TestFixture::NodeT;

    const auto G = this->unweightedGraph();
    std::vector<NodeT> sequence;
    Traversal::DFSfrom(G, NodeT{0}, [&](auto u) {
        static_assert(std::is_same_v<decltype(u), NodeT>);
        sequence.push_back(u);
    });

    EXPECT_THAT(sequence,
                ElementsAre(NodeT{0}, NodeT{2}, NodeT{4}, NodeT{6}, NodeT{1}, NodeT{3}, NodeT{5}));
}

TYPED_TEST(GenericDFSTraversalGTest, testDFSfromWeightedGraph) {
    using NodeT = typename TestFixture::NodeT;

    const auto G = this->weightedGraph();
    std::vector<NodeT> sequence;
    Traversal::DFSfrom(G, NodeT{0}, [&](NodeT u) { sequence.push_back(u); });

    EXPECT_THAT(sequence,
                ElementsAre(NodeT{0}, NodeT{2}, NodeT{3}, NodeT{5}, NodeT{4}, NodeT{6}, NodeT{1}));
}

TYPED_TEST(GenericDFSTraversalGTest, testDFSfromDirectedGraph) {
    using NodeT = typename TestFixture::NodeT;

    const auto G = this->directedGraph();
    std::vector<NodeT> sequence;
    Traversal::DFSfrom(G, NodeT{0}, [&](NodeT u) { sequence.push_back(u); });

    EXPECT_THAT(sequence, ElementsAre(NodeT{0}, NodeT{1}, NodeT{2}, NodeT{4}));
}

TYPED_TEST(GenericDFSTraversalGTest, testDFSfromIsolatedNode) {
    using NodeT = typename TestFixture::NodeT;

    const auto G = this->directedGraph();
    std::vector<NodeT> sequence;
    Traversal::DFSfrom(G, NodeT{5}, [&](NodeT u) { sequence.push_back(u); });

    EXPECT_THAT(sequence, ElementsAre(NodeT{5}));
}

TYPED_TEST(GenericDFSTraversalGTest, testDFSfromRemovedNodes) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    auto G = this->pathGraph();
    G.removeNode(NodeT{1});

    std::vector<NodeT> sequence;
    Traversal::DFSfrom(G, NodeT{4}, [&](NodeT u) { sequence.push_back(u); });
    EXPECT_THAT(sequence, ElementsAre(NodeT{4}, NodeT{3}, NodeT{2}));

    std::vector<std::pair<NodeT, NodeT>> edgeSequence;
    Traversal::DFSEdgesFrom(G, NodeT{4}, [&](NodeT u, NodeT v, EdgeWeightT, edgeid) {
        edgeSequence.emplace_back(u, v);
    });
    EXPECT_THAT(edgeSequence, ElementsAre(Pair(NodeT{4}, NodeT{3}), Pair(NodeT{3}, NodeT{2})));
}

TYPED_TEST(GenericDFSTraversalGTest, testDFSfromInvalidSource) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;

    auto G = this->pathGraph();
    count calls = 0;
    auto handle = [&](NodeT) { ++calls; };
    EXPECT_THROW(Traversal::DFSfrom(G, static_cast<NodeT>(G.upperNodeIdBound()), handle),
                 std::runtime_error);

    G.removeNode(NodeT{2});
    EXPECT_THROW(Traversal::DFSfrom(G, NodeT{2}, handle), std::runtime_error);
    if constexpr (std::is_signed_v<NodeT>) {
        EXPECT_THROW(Traversal::DFSfrom(G, NodeT{-1}, handle), std::runtime_error);
    }
    EXPECT_THROW(Traversal::DFSfrom(GraphT(0), NodeT{0}, handle), std::runtime_error);
    EXPECT_EQ(calls, 0);
}

TYPED_TEST(GenericDFSTraversalGTest, testDFSEdgesFromTypedGraphs) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    auto G = this->weightedGraph();
    G.indexEdges();
    std::vector<std::tuple<NodeT, NodeT, EdgeWeightT, edgeid>> edgeSequence;

    Traversal::DFSEdgesFrom(G, NodeT{0}, [&](auto u, auto v, auto w, auto eid) {
        static_assert(std::is_same_v<decltype(u), NodeT>);
        static_assert(std::is_same_v<decltype(v), NodeT>);
        static_assert(std::is_same_v<decltype(w), EdgeWeightT>);
        static_assert(std::is_same_v<decltype(eid), edgeid>);
        edgeSequence.emplace_back(u, v, w, eid);
    });

    EXPECT_THAT(
        edgeSequence,
        ElementsAre(FieldsAre(NodeT{0}, NodeT{1}, EdgeWeightT{10}, G.edgeId(NodeT{0}, NodeT{1})),
                    FieldsAre(NodeT{0}, NodeT{2}, EdgeWeightT{20}, G.edgeId(NodeT{0}, NodeT{2})),
                    FieldsAre(NodeT{2}, NodeT{4}, EdgeWeightT{40}, G.edgeId(NodeT{2}, NodeT{4})),
                    FieldsAre(NodeT{2}, NodeT{3}, EdgeWeightT{70}, G.edgeId(NodeT{2}, NodeT{3})),
                    FieldsAre(NodeT{3}, NodeT{5}, EdgeWeightT{50}, G.edgeId(NodeT{3}, NodeT{5})),
                    FieldsAre(NodeT{4}, NodeT{6}, EdgeWeightT{60}, G.edgeId(NodeT{4}, NodeT{6}))));
}

TYPED_TEST(GenericDFSTraversalGTest, testDFSEdgesFromCornerCases) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    GraphT G(5, true, true);
    G.addEdge(NodeT{0}, NodeT{0}, EdgeWeightT{10});
    G.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{20});
    G.addEdge(NodeT{0}, NodeT{2}, EdgeWeightT{30});
    G.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{40});
    G.addEdge(NodeT{2}, NodeT{0}, EdgeWeightT{50});
    G.addEdge(NodeT{3}, NodeT{0}, EdgeWeightT{60});
    G.indexEdges();

    std::vector<std::tuple<NodeT, NodeT, EdgeWeightT, edgeid>> edgeSequence;
    auto handle = [&](NodeT u, NodeT v, EdgeWeightT w, edgeid eid) {
        edgeSequence.emplace_back(u, v, w, eid);
    };
    Traversal::DFSEdgesFrom(G, NodeT{0}, handle);

    EXPECT_THAT(
        edgeSequence,
        ElementsAre(FieldsAre(NodeT{0}, NodeT{1}, EdgeWeightT{20}, G.edgeId(NodeT{0}, NodeT{1})),
                    FieldsAre(NodeT{0}, NodeT{2}, EdgeWeightT{30}, G.edgeId(NodeT{0}, NodeT{2}))));

    edgeSequence.clear();
    Traversal::DFSEdgesFrom(G, NodeT{4}, handle);
    EXPECT_THAT(edgeSequence, IsEmpty());
}

TYPED_TEST(GenericDFSTraversalGTest, testDFSEdgesFromInvalidSource) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    auto G = this->weightedGraph();
    count calls = 0;
    auto handle = [&](NodeT, NodeT, EdgeWeightT, edgeid) { ++calls; };
    EXPECT_THROW(Traversal::DFSEdgesFrom(G, static_cast<NodeT>(G.upperNodeIdBound()), handle),
                 std::runtime_error);

    G.removeNode(NodeT{2});
    EXPECT_THROW(Traversal::DFSEdgesFrom(G, NodeT{2}, handle), std::runtime_error);
    if constexpr (std::is_signed_v<NodeT>) {
        EXPECT_THROW(Traversal::DFSEdgesFrom(G, NodeT{-1}, handle), std::runtime_error);
    }
    EXPECT_THROW(Traversal::DFSEdgesFrom(GraphT(0), NodeT{0}, handle), std::runtime_error);
    EXPECT_EQ(calls, 0);
}

} // namespace
} // namespace NetworKit
