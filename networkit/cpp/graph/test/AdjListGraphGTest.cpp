#include <gtest/gtest.h>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace {

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

REGISTER_TYPED_TEST_SUITE_P(AdjListGraphGTest, testDefaultConstructor);

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
