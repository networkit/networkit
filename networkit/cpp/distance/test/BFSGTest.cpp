#include <span>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/distance/BFS.hpp>
#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/graph/GraphTestUtils.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {
namespace {

using ::testing::AllOf;
using ::testing::Each;
using ::testing::Eq;
using ::testing::HasSubstr;
using ::testing::IsEmpty;
using ::testing::Not;
using ::testing::ResultOf;
using ::testing::ThrowsMessage;

template <class NodeT_, class EdgeWeightT_>
struct BFSConfig {
    using NodeT = NodeT_;
    using EdgeWeightT = EdgeWeightT_;
};

template <class TestT>
class BFSGTest : public ::testing::Test {
public:
    using NodeT = typename TestT::NodeT;
    using EdgeWeightT = typename TestT::EdgeWeightT;
};

TYPED_TEST_SUITE_P(BFSGTest);

TYPED_TEST_P(BFSGTest, testBFSThrowsInvalidSource) {
    using GraphT = AdjListGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT>;
    using NodeT = typename TestFixture::NodeT;
    GraphT G(5, true);
    BreadthFirstSearch<GraphT> bfs(G, NodeT{1});
    EXPECT_THAT([&]() { bfs.setSource(NodeT{6}); },
                ThrowsMessage<std::runtime_error>(HasSubstr("node not in the graph")));
}

TYPED_TEST_P(BFSGTest, testBFSThrowsInvalidTarget) {
    using GraphT = AdjListGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT>;
    using NodeT = typename TestFixture::NodeT;
    GraphT G(5, true);
    BreadthFirstSearch<GraphT> bfs(G, NodeT{1});
    EXPECT_THAT([&]() { bfs.setTarget(NodeT{6}); },
                ThrowsMessage<std::runtime_error>(HasSubstr("node not in the graph")));
}

TYPED_TEST_P(BFSGTest, testBFSThrowsStorepathNotTrueCallPredecessor) {
    using GraphT = AdjListGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT>;
    using NodeT = typename TestFixture::NodeT;
    GraphT G(3, true);
    G.addEdge(NodeT{0}, NodeT{2});
    G.addEdge(NodeT{1}, NodeT{2});
    BreadthFirstSearch<GraphT> bfs(G, NodeT{0}, /*storePaths*/ false);
    bfs.setTarget(NodeT{2});
    bfs.run();
    EXPECT_THAT([&]() { bfs.getPredecessors(NodeT{1}); },
                ThrowsMessage<std::runtime_error>(HasSubstr("predecessors have not been stored")));
}

TYPED_TEST_P(BFSGTest, testShortestPaths) {
    using GraphT = AdjListGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT>;
    using NodeT = typename TestFixture::NodeT;
    const GraphT G_typed =
        TestUtils::ToTGraph<GraphT>(METISGraphReader{}.read("input/PGPgiantcompo.graph"));
    const NodeT source{2};
    BreadthFirstSearch<GraphT> bfs(G_typed, source);
    bfs.run();
    bigfloat max = 0;
    NodeT x{0};
    G_typed.forNodes([&](NodeT n) -> void {
        if (bfs.numberOfPaths(n) > max) {
            max = bfs.numberOfPaths(n);
            x = n;
        }
    });
    const count dist = bfs.distance(x);
    EXPECT_THAT(
        bfs.getPaths(x, true),
        AllOf(Not(IsEmpty()),
              Each(AllOf(ResultOf([](std::span<const NodeT> p) { return p[0]; }, Eq(source)),
                         ResultOf([dist](std::span<const NodeT> p) { return p[dist]; }, Eq(x))))));
}

TYPED_TEST_P(BFSGTest, testDirectedBFS) {
    using GraphT = AdjListGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT>;
    using NodeT = typename TestFixture::NodeT;
    GraphT G(8, false, true);

    G.addEdge(NodeT{0}, NodeT{6});
    G.addEdge(NodeT{0}, NodeT{2});
    G.addEdge(NodeT{3}, NodeT{2});
    G.addEdge(NodeT{5}, NodeT{3});
    G.addEdge(NodeT{6}, NodeT{5});
    G.addEdge(NodeT{5}, NodeT{7});
    G.addEdge(NodeT{4}, NodeT{5});
    G.addEdge(NodeT{2}, NodeT{4});
    G.addEdge(NodeT{2}, NodeT{1});

    BreadthFirstSearch<GraphT> sssp(G, NodeT{0});
    sssp.run();
    EXPECT_EQ(sssp.distance(NodeT{0}), 0);
    EXPECT_EQ(sssp.distance(NodeT{1}), 2);
    EXPECT_EQ(sssp.distance(NodeT{2}), 1);
    EXPECT_EQ(sssp.distance(NodeT{3}), 3);
    EXPECT_EQ(sssp.distance(NodeT{4}), 2);
    EXPECT_EQ(sssp.distance(NodeT{5}), 2);
    EXPECT_EQ(sssp.distance(NodeT{6}), 1);
    EXPECT_EQ(sssp.distance(NodeT{7}), 3);
}

REGISTER_TYPED_TEST_SUITE_P(BFSGTest, testBFSThrowsInvalidSource, testBFSThrowsInvalidTarget,
                            testBFSThrowsStorepathNotTrueCallPredecessor, testShortestPaths,
                            testDirectedBFS);

using BFSTestTypes =
    ::testing::Types<BFSConfig<node, edgeweight>, BFSConfig<int, float>, BFSConfig<int, int>>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestBFS, BFSGTest, BFSTestTypes, );

} // namespace
} // namespace NetworKit
