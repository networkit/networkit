#include <gtest/gtest.h>
#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/graph/GraphTestUtils.hpp>

namespace NetworKit::TestUtils {

TEST(GraphTestUtils, SimpleConversion) {
    AdjListGraph<node, edgeweight> G(3, true, false);
    G.addEdge(0, 1, 1.0);
    G.addEdge(1, 2, 2.5);

    auto G_typed = ToTGraph<AdjListGraph<uint32_t, double>>(G);
    EXPECT_EQ(G.numberOfNodes(), G_typed.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), G_typed.numberOfEdges());
    EXPECT_TRUE(G_typed.hasEdge(0, 1));
    EXPECT_TRUE(G_typed.hasEdge(1, 2));
    EXPECT_DOUBLE_EQ(G_typed.weight(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(G_typed.weight(1, 2), 2.5);
}

} // namespace NetworKit::TestUtils
