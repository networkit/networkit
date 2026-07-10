#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace {

MATCHER_P2(EdgeEq, expectedU, expectedV, "") {
    return arg.u == expectedU && arg.v == expectedV;
}

MATCHER_P3(WeightedEdgeEq, expectedU, expectedV, expectedWeight, "") {
    return arg.u == expectedU && arg.v == expectedV && arg.weight == expectedWeight;
}

/**
 * Configuration for AdjListGraph tests, combining type and value parameterization.
 */
template <class NodeT_, class EdgeWeightT_, bool Weighted, bool Directed>
struct AdjListConfig {
    using NodeT = NodeT_;
    using EdgeWeightT = EdgeWeightT_;
    static constexpr bool weighted = Weighted;
    static constexpr bool directed = Directed;
};

template <class TestT>
class AdjListGraphGTest : public ::testing::Test {
public:
    using NodeT = typename TestT::NodeT;
    using EdgeWeightT = typename TestT::EdgeWeightT;
};

TYPED_TEST_SUITE_P(AdjListGraphGTest);

TYPED_TEST_P(AdjListGraphGTest, testDefaultConstructor) {
    AdjListGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT> G(
        10, TypeParam::weighted, TypeParam::directed);
    EXPECT_EQ(G.numberOfNodes(), 10u);
    EXPECT_EQ(G.numberOfEdges(), 0u);
    EXPECT_EQ(G.isWeighted(), TypeParam::weighted);
    EXPECT_EQ(G.isDirected(), TypeParam::directed);
}

TYPED_TEST_P(AdjListGraphGTest, testNodeAndEdgeIterators) {
    using ::testing::ElementsAre;

    using Graph = AdjListGraph<typename TestFixture::NodeT, typename TestFixture::EdgeWeightT>;
    using NodeT = typename TestFixture::NodeT;
    using EdgeWeightT = typename TestFixture::EdgeWeightT;

    Graph G(4, TypeParam::weighted, TypeParam::directed);
    G.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{2});
    G.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{3});

    const std::vector<NodeT> nodes(G.nodeRange().begin(), G.nodeRange().end());
    EXPECT_THAT(nodes, ElementsAre(NodeT{0}, NodeT{1}, NodeT{2}, NodeT{3}));

    const std::vector<EdgeT<NodeT>> edges(G.edgeRange().begin(), G.edgeRange().end());
    EXPECT_THAT(edges, ElementsAre(EdgeEq(NodeT{0}, NodeT{1}), EdgeEq(NodeT{1}, NodeT{2})));

    const EdgeWeightT firstWeight = TypeParam::weighted ? EdgeWeightT{2} : EdgeWeightT{1};
    const EdgeWeightT secondWeight = TypeParam::weighted ? EdgeWeightT{3} : EdgeWeightT{1};
    const std::vector<WeightedEdgeT<NodeT, EdgeWeightT>> weightedEdges(G.edgeWeightRange().begin(),
                                                                       G.edgeWeightRange().end());
    EXPECT_THAT(weightedEdges, ElementsAre(WeightedEdgeEq(NodeT{0}, NodeT{1}, firstWeight),
                                           WeightedEdgeEq(NodeT{1}, NodeT{2}, secondWeight)));
}

REGISTER_TYPED_TEST_SUITE_P(AdjListGraphGTest, testDefaultConstructor, testNodeAndEdgeIterators);

using AdjListTestTypes = ::testing::Types<
    AdjListConfig<node, edgeweight, false, false>, AdjListConfig<node, edgeweight, true, false>,
    AdjListConfig<node, edgeweight, false, true>, AdjListConfig<node, edgeweight, true, true>,
    AdjListConfig<int, float, false, false>, AdjListConfig<int, float, true, false>,
    AdjListConfig<int, float, false, true>, AdjListConfig<int, float, true, true>,
    AdjListConfig<int, int, false, false>, AdjListConfig<int, int, true, false>,
    AdjListConfig<int, int, false, true>, AdjListConfig<int, int, true, true>>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestAdjListGraph, AdjListGraphGTest, AdjListTestTypes, );

} // namespace
} // namespace NetworKit
