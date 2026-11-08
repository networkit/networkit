/* DijkstraTraversalGTest.cpp */

#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/Dijkstra.hpp>

#include "TraversalTestHelper.hpp"

namespace NetworKit {
namespace {

using ::testing::ElementsAre;
using ::testing::IsEmpty;
using ::testing::Pair;
using ::testing::UnorderedElementsAre;

template <class TestT>
using GenericDijkstraTraversalGTest = GenericTraversalGTest<TestT>;

// Includes graphs with identical node and weight types.
TYPED_TEST_SUITE(GenericDijkstraTraversalGTest, TraversalTestTypes,
                 /*Comma needed for variadic macro.*/);

TYPED_TEST(GenericDijkstraTraversalGTest, testDijkstraFromTypedGraphs) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    const auto G = this->unweightedGraph();
    std::vector<std::pair<NodeT, EdgeWeightT>> visits;
    std::vector<EdgeWeightT> distances;
    Traversal::DijkstraFrom(G, NodeT{0}, [&](auto u, auto dist) {
        static_assert(std::is_same_v<decltype(u), NodeT>);
        static_assert(std::is_same_v<decltype(dist), EdgeWeightT>);
        visits.emplace_back(u, dist);
        distances.push_back(dist);
    });

    // Equal-distance nodes may be visited in either order.
    EXPECT_THAT(visits,
                UnorderedElementsAre(Pair(NodeT{0}, EdgeWeightT{0}), Pair(NodeT{1}, EdgeWeightT{1}),
                                     Pair(NodeT{2}, EdgeWeightT{1}), Pair(NodeT{3}, EdgeWeightT{2}),
                                     Pair(NodeT{4}, EdgeWeightT{2}), Pair(NodeT{5}, EdgeWeightT{3}),
                                     Pair(NodeT{6}, EdgeWeightT{3})));
    EXPECT_TRUE(std::is_sorted(distances.begin(), distances.end()));
}

TYPED_TEST(GenericDijkstraTraversalGTest, testDijkstraFromWeightedGraph) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    const auto G = this->weightedGraph();
    std::vector<NodeT> sequence;
    std::vector<EdgeWeightT> distances;
    Traversal::DijkstraFrom(G, NodeT{0}, [&](NodeT u, EdgeWeightT dist) {
        sequence.push_back(u);
        distances.push_back(dist);
    });

    EXPECT_THAT(sequence,
                ElementsAre(NodeT{0}, NodeT{1}, NodeT{2}, NodeT{3}, NodeT{4}, NodeT{5}, NodeT{6}));
    EXPECT_THAT(distances, ElementsAre(0, 10, 20, 40, 60, 90, 120));
}

TYPED_TEST(GenericDijkstraTraversalGTest, testDijkstraFromDistanceUpdates) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    GraphT G(4, true, true);
    G.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{10});
    G.addEdge(NodeT{0}, NodeT{2}, EdgeWeightT{1});
    G.addEdge(NodeT{2}, NodeT{1}, EdgeWeightT{1});
    G.addEdge(NodeT{1}, NodeT{3}, EdgeWeightT{1});
    G.addEdge(NodeT{2}, NodeT{3}, EdgeWeightT{20});

    std::vector<std::pair<NodeT, EdgeWeightT>> visits;
    Traversal::DijkstraFrom(G, NodeT{0},
                            [&](NodeT u, EdgeWeightT dist) { visits.emplace_back(u, dist); });

    // Improved nodes must be reordered and visited exactly once, including when
    // the queue implementation retains outdated entries.
    EXPECT_THAT(visits,
                ElementsAre(Pair(NodeT{0}, EdgeWeightT{0}), Pair(NodeT{2}, EdgeWeightT{1}),
                            Pair(NodeT{1}, EdgeWeightT{2}), Pair(NodeT{3}, EdgeWeightT{3})));
}

TYPED_TEST(GenericDijkstraTraversalGTest, testDijkstraFromMultiSourceDistances) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    const auto G = this->pathGraph();
    const std::vector<NodeT> sources{NodeT{0}, NodeT{4}};
    std::vector<std::pair<NodeT, EdgeWeightT>> visits;
    std::vector<EdgeWeightT> distances;
    Traversal::DijkstraFrom(G, sources.begin(), sources.end(), [&](NodeT u, EdgeWeightT dist) {
        visits.emplace_back(u, dist);
        distances.push_back(dist);
    });

    EXPECT_THAT(visits,
                UnorderedElementsAre(Pair(NodeT{0}, EdgeWeightT{0}), Pair(NodeT{4}, EdgeWeightT{0}),
                                     Pair(NodeT{1}, EdgeWeightT{1}), Pair(NodeT{3}, EdgeWeightT{1}),
                                     Pair(NodeT{2}, EdgeWeightT{2})));
    EXPECT_TRUE(std::is_sorted(distances.begin(), distances.end()));
}

TYPED_TEST(GenericDijkstraTraversalGTest, testDijkstraFromEmptySourceRange) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;

    const std::vector<NodeT> sources;
    std::vector<NodeT> sequence;
    auto handle = [&](NodeT u) { sequence.push_back(u); };
    Traversal::DijkstraFrom(this->pathGraph(), sources.begin(), sources.end(), handle);
    Traversal::DijkstraFrom(GraphT(0), sources.begin(), sources.end(), handle);
    EXPECT_THAT(sequence, IsEmpty());
}

TYPED_TEST(GenericDijkstraTraversalGTest, testDijkstraFromDuplicateSources) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    const auto G = this->pathGraph();
    const std::vector<NodeT> sources{NodeT{0}, NodeT{0}, NodeT{2}};
    std::vector<std::pair<NodeT, EdgeWeightT>> visits;
    Traversal::DijkstraFrom(G, sources.begin(), sources.end(),
                            [&](NodeT u, EdgeWeightT dist) { visits.emplace_back(u, dist); });

    EXPECT_THAT(visits,
                UnorderedElementsAre(Pair(NodeT{0}, EdgeWeightT{0}), Pair(NodeT{2}, EdgeWeightT{0}),
                                     Pair(NodeT{1}, EdgeWeightT{1}), Pair(NodeT{3}, EdgeWeightT{1}),
                                     Pair(NodeT{4}, EdgeWeightT{2})));
}

TYPED_TEST(GenericDijkstraTraversalGTest, testDijkstraFromDirectedGraph) {
    using NodeT = typename TestFixture::NodeT;

    const auto G = this->directedGraph();
    std::vector<NodeT> sequence;
    Traversal::DijkstraFrom(G, NodeT{0}, [&](auto u) {
        static_assert(std::is_same_v<decltype(u), NodeT>);
        sequence.push_back(u);
    });
    EXPECT_THAT(sequence, ElementsAre(NodeT{0}, NodeT{1}, NodeT{2}, NodeT{4}));

    sequence.clear();
    Traversal::DijkstraFrom(G, NodeT{5}, [&](NodeT u) { sequence.push_back(u); });
    EXPECT_THAT(sequence, ElementsAre(NodeT{5}));
}

TYPED_TEST(GenericDijkstraTraversalGTest, testDijkstraFromRemovedNodes) {
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    auto G = this->pathGraph();
    G.removeNode(NodeT{1});
    std::vector<std::pair<NodeT, EdgeWeightT>> visits;
    Traversal::DijkstraFrom(G, NodeT{4},
                            [&](NodeT u, EdgeWeightT dist) { visits.emplace_back(u, dist); });
    EXPECT_THAT(visits, ElementsAre(Pair(NodeT{4}, EdgeWeightT{0}), Pair(NodeT{3}, EdgeWeightT{1}),
                                    Pair(NodeT{2}, EdgeWeightT{2})));
}

TYPED_TEST(GenericDijkstraTraversalGTest, testDijkstraFromZeroWeightCycle) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    GraphT G(4, true, true);
    G.addEdge(NodeT{0}, NodeT{0}, EdgeWeightT{0});
    G.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{0});
    G.addEdge(NodeT{1}, NodeT{0}, EdgeWeightT{0});
    G.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{2});
    std::vector<std::pair<NodeT, EdgeWeightT>> visits;
    Traversal::DijkstraFrom(G, NodeT{0},
                            [&](NodeT u, EdgeWeightT dist) { visits.emplace_back(u, dist); });
    EXPECT_THAT(visits, ElementsAre(Pair(NodeT{0}, EdgeWeightT{0}), Pair(NodeT{1}, EdgeWeightT{0}),
                                    Pair(NodeT{2}, EdgeWeightT{2})));
}

TYPED_TEST(GenericDijkstraTraversalGTest, testDijkstraFromInvalidSource) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;

    auto G = this->pathGraph();
    count calls = 0;
    auto handle = [&](NodeT) { ++calls; };
    EXPECT_THROW(Traversal::DijkstraFrom(G, static_cast<NodeT>(G.upperNodeIdBound()), handle),
                 std::runtime_error);
    G.removeNode(NodeT{2});
    EXPECT_THROW(Traversal::DijkstraFrom(G, NodeT{2}, handle), std::runtime_error);
    if constexpr (std::is_signed_v<NodeT>) {
        EXPECT_THROW(Traversal::DijkstraFrom(G, NodeT{-1}, handle), std::runtime_error);
    }
    EXPECT_THROW(Traversal::DijkstraFrom(GraphT(0), NodeT{0}, handle), std::runtime_error);
    EXPECT_EQ(calls, 0);
}

TYPED_TEST(GenericDijkstraTraversalGTest, testDijkstraFromInvalidSources) {
    using NodeT = typename TestFixture::NodeT;

    auto G = this->pathGraph();
    count calls = 0;
    auto handle = [&](NodeT) { ++calls; };
    const std::vector<NodeT> outOfRange{NodeT{0}, static_cast<NodeT>(G.upperNodeIdBound())};
    EXPECT_THROW(Traversal::DijkstraFrom(G, outOfRange.begin(), outOfRange.end(), handle),
                 std::runtime_error);
    G.removeNode(NodeT{2});
    const std::vector<NodeT> removed{NodeT{0}, NodeT{2}};
    EXPECT_THROW(Traversal::DijkstraFrom(G, removed.begin(), removed.end(), handle),
                 std::runtime_error);
    if constexpr (std::is_signed_v<NodeT>) {
        const std::vector<NodeT> negative{NodeT{0}, NodeT{-1}};
        EXPECT_THROW(Traversal::DijkstraFrom(G, negative.begin(), negative.end(), handle),
                     std::runtime_error);
    }
    EXPECT_EQ(calls, 0);
}

} // namespace
} // namespace NetworKit
