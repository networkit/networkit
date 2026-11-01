#include <array>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/graph/SteinerTree.hpp>
#include <networkit/structures/UnionFind.hpp>

namespace NetworKit {

namespace {

template <class NodeT_, class EdgeWeightT_>
struct SteinerTreeConfig {
    using NodeT = NodeT_;
    using EdgeWeightT = EdgeWeightT_;
};

template <class TestT>
class SteinerTreeGTest : public testing::Test {
public:
    using NodeT = typename TestT::NodeT;
    using EdgeWeightT = typename TestT::EdgeWeightT;
    using GraphT = AdjListGraph<NodeT, EdgeWeightT>;
    using SteinerTreeT = SteinerTree<GraphT>;

    static void expectValidTree(const GraphT &G, const GraphT &tree,
                                std::span<const NodeT> terminals) {
        UnionFind components(G.upperNodeIdBound());
        uint64_t nodesInTree = 0;
        tree.forNodes([&](NodeT u) -> void {
            if (tree.degree(u) > 0) {
                ++nodesInTree;
            }
            // Non-terminal leaves would only add weight.
            if (tree.degree(u) == 1) {
                EXPECT_THAT(terminals, testing::Contains(u));
            }
        });
        tree.forEdges([&](NodeT u, NodeT v, EdgeWeightT weight, edgeid) -> void {
            EXPECT_TRUE(G.hasEdge(u, v));
            EXPECT_EQ(G.weight(u, v), weight);
            EXPECT_NE(components.find(u), components.find(v)) << "cycle through " << u;
            components.merge(u, v);
        });
        if (terminals.size() > 1) {
            EXPECT_EQ(tree.numberOfEdges() + 1, nodesInTree);
        }
        for (const NodeT terminal : terminals) {
            EXPECT_EQ(components.find(terminal), components.find(terminals.front()));
        }
    }
};

TYPED_TEST_SUITE_P(SteinerTreeGTest);

TYPED_TEST_P(SteinerTreeGTest, testUsesASteinerNode) {
    using SteinerTreeT = typename TestFixture::SteinerTreeT;

    // A star is lighter than any tree on the terminals alone.
    typename TestFixture::GraphT G(4, true);
    G.addEdge(0, 3, 1);
    G.addEdge(1, 3, 1);
    G.addEdge(2, 3, 1);
    G.addEdge(0, 1, 3);
    G.addEdge(1, 2, 3);
    const auto terminals = std::to_array<typename TestFixture::NodeT>({0, 1, 2});

    SteinerTreeT algo(G, terminals);
    algo.run();

    ASSERT_TRUE(algo.hasTree());
    EXPECT_EQ(algo.getTotalWeight(), typename SteinerTreeT::WeightSumT{3});
    EXPECT_EQ(algo.getTree().degree(3), 3);
    TestFixture::expectValidTree(G, algo.getTree(), terminals);
}

TYPED_TEST_P(SteinerTreeGTest, testPathBetweenTwoTerminalsIsShortest) {
    using SteinerTreeT = typename TestFixture::SteinerTreeT;

    typename TestFixture::GraphT G(5, true);
    G.addEdge(0, 1, 1);
    G.addEdge(1, 2, 1);
    G.addEdge(2, 3, 1);
    G.addEdge(0, 4, 2);
    G.addEdge(4, 3, 5);
    const auto terminals = std::to_array<typename TestFixture::NodeT>({0, 3});

    SteinerTreeT algo(G, terminals);
    algo.run();

    EXPECT_EQ(algo.getTotalWeight(), typename SteinerTreeT::WeightSumT{3});
    TestFixture::expectValidTree(G, algo.getTree(), terminals);
}

TYPED_TEST_P(SteinerTreeGTest, testInducedTreeIsImprovedAndPruned) {
    using SteinerTreeT = typename TestFixture::SteinerTreeT;

    // Nodes 3 and 6 are closest to terminal 0, node 4 to terminal 1 and node 5 to terminal 2. The
    // shortest bridges are 5-0 (length 34) and 6-4 (length 42), so Mehlhorn's tree is
    // 2-5-0-3-6-4-1 with weight 76. The minimum spanning tree of all seven nodes replaces 6-4 by
    // 5-4 (weight 74), which leaves the non-terminals 3 and 6 hanging from 0. Removing 6 makes 3 a
    // leaf, and removing it too gives the optimal tree.
    typename TestFixture::GraphT G(7, true);
    G.addEdge(0, 3, 5);
    G.addEdge(3, 6, 5);
    G.addEdge(6, 4, 20);
    G.addEdge(4, 1, 12);
    G.addEdge(2, 5, 15);
    G.addEdge(5, 4, 18);
    G.addEdge(5, 0, 19);
    const auto terminals = std::to_array<typename TestFixture::NodeT>({0, 1, 2});

    SteinerTreeT algo(G, terminals);
    algo.run();

    ASSERT_TRUE(algo.hasTree());
    EXPECT_EQ(algo.getTotalWeight(), typename SteinerTreeT::WeightSumT{64});
    EXPECT_EQ(algo.getTree().numberOfEdges(), 4);
    EXPECT_EQ(algo.getTree().degree(3), 0);
    EXPECT_EQ(algo.getTree().degree(6), 0);
    EXPECT_TRUE(algo.getTree().hasEdge(5, 4));
    TestFixture::expectValidTree(G, algo.getTree(), terminals);
}

TYPED_TEST_P(SteinerTreeGTest, testLongPathsDoNotOverflow) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using SteinerTreeT = typename TestFixture::SteinerTreeT;
    using WeightSumT = typename SteinerTreeT::WeightSumT;

    // Every weight fits in uint8_t, but path lengths do not. The bridge 0-1 has length 300, which
    // must not be taken for shorter than the bridge on the path 0-2-3 of length 100.
    GraphT G(4, true);
    G.addEdge(0, 1, 200);
    G.addEdge(1, 3, 100);
    G.addEdge(0, 2, 50);
    G.addEdge(2, 3, 50);
    const auto terminals = std::to_array<NodeT>({0, 3});

    SteinerTreeT algo(G, terminals);
    algo.run();

    EXPECT_EQ(algo.getTotalWeight(), WeightSumT{100});
    EXPECT_TRUE(algo.getTree().hasEdge(0, 2));
    TestFixture::expectValidTree(G, algo.getTree(), terminals);

    // Node 1 is at distance 255 from both terminals, the largest uint8_t, and the tree weighs 510.
    GraphT path(3, true);
    path.addEdge(0, 1, 255);
    path.addEdge(1, 2, 255);
    const auto pathTerminals = std::to_array<NodeT>({0, 2});

    SteinerTreeT pathAlgo(path, pathTerminals);
    pathAlgo.run();

    ASSERT_TRUE(pathAlgo.hasTree());
    EXPECT_EQ(pathAlgo.getTotalWeight(), WeightSumT{510});
    TestFixture::expectValidTree(path, pathAlgo.getTree(), pathTerminals);
}

TYPED_TEST_P(SteinerTreeGTest, testUnweightedGraph) {
    using SteinerTreeT = typename TestFixture::SteinerTreeT;

    typename TestFixture::GraphT G(5);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(1, 3);
    G.addEdge(3, 4);
    const auto terminals = std::to_array<typename TestFixture::NodeT>({0, 2, 4});

    SteinerTreeT algo(G, terminals);
    algo.run();

    EXPECT_EQ(algo.getTotalWeight(), typename SteinerTreeT::WeightSumT{4});
    EXPECT_TRUE(algo.getTree().isWeighted());
    TestFixture::expectValidTree(G, algo.getTree(), terminals);
}

TYPED_TEST_P(SteinerTreeGTest, testAllNodesAreTerminals) {
    using SteinerTreeT = typename TestFixture::SteinerTreeT;

    // Every node is its own region, so the bridges are the edges and the result is a minimum
    // spanning tree. The bridge 2-0 is skipped because 0 and 2 are already connected through 1.
    typename TestFixture::GraphT G(4, true);
    G.addEdge(0, 1, 1);
    G.addEdge(1, 2, 2);
    G.addEdge(2, 0, 3);
    G.addEdge(2, 3, 4);
    const auto terminals = std::to_array<typename TestFixture::NodeT>({0, 1, 2, 3});

    SteinerTreeT algo(G, terminals);
    algo.run();

    EXPECT_EQ(algo.getTotalWeight(), typename SteinerTreeT::WeightSumT{7});
    EXPECT_FALSE(algo.getTree().hasEdge(2, 0));
    TestFixture::expectValidTree(G, algo.getTree(), terminals);
}

TYPED_TEST_P(SteinerTreeGTest, testDisconnectedTerminals) {
    using SteinerTreeT = typename TestFixture::SteinerTreeT;

    typename TestFixture::GraphT G(4);
    G.addEdge(0, 1);
    G.addEdge(2, 3);

    SteinerTreeT algo(G, std::to_array<typename TestFixture::NodeT>({0, 1, 2}));
    algo.run();

    EXPECT_FALSE(algo.hasTree());
    EXPECT_EQ(algo.getTree().numberOfEdges(), 0);
    EXPECT_EQ(algo.getTotalWeight(), typename SteinerTreeT::WeightSumT{0});
}

TYPED_TEST_P(SteinerTreeGTest, testComponentsWithoutTerminalsAreIgnored) {
    using SteinerTreeT = typename TestFixture::SteinerTreeT;

    typename TestFixture::GraphT G(5);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(3, 4);

    SteinerTreeT algo(G, std::to_array<typename TestFixture::NodeT>({0, 2}));
    algo.run();

    ASSERT_TRUE(algo.hasTree());
    EXPECT_EQ(algo.getTotalWeight(), typename SteinerTreeT::WeightSumT{2});
}

TYPED_TEST_P(SteinerTreeGTest, testFewTerminals) {
    using NodeT = typename TestFixture::NodeT;
    using SteinerTreeT = typename TestFixture::SteinerTreeT;

    typename TestFixture::GraphT G(3);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    for (const std::vector<NodeT> &terminals :
         {std::vector<NodeT>{}, std::vector<NodeT>{1}, std::vector<NodeT>{1, 1}}) {
        SteinerTreeT algo(G, terminals);
        algo.run();
        ASSERT_TRUE(algo.hasTree());
        EXPECT_EQ(algo.getTree().numberOfEdges(), 0);
        EXPECT_EQ(algo.getTotalWeight(), typename SteinerTreeT::WeightSumT{0});
    }
}

TYPED_TEST_P(SteinerTreeGTest, testZeroWeightsAndSelfLoops) {
    using SteinerTreeT = typename TestFixture::SteinerTreeT;

    typename TestFixture::GraphT G(4, true);
    G.addEdge(0, 1, 0);
    G.addEdge(1, 2, 0);
    G.addEdge(1, 1, 3);
    G.addEdge(2, 3, 1);
    const auto terminals = std::to_array<typename TestFixture::NodeT>({0, 3});

    SteinerTreeT algo(G, terminals);
    algo.run();

    EXPECT_EQ(algo.getTotalWeight(), typename SteinerTreeT::WeightSumT{1});
    TestFixture::expectValidTree(G, algo.getTree(), terminals);
}

TYPED_TEST_P(SteinerTreeGTest, testDeletedNodes) {
    using SteinerTreeT = typename TestFixture::SteinerTreeT;

    typename TestFixture::GraphT G(5, true);
    G.addEdge(0, 2, 1);
    G.addEdge(2, 4, 1);
    G.addEdge(1, 0, 1);
    G.addEdge(1, 4, 1);
    G.removeNode(1);
    G.removeNode(3);
    const auto terminals = std::to_array<typename TestFixture::NodeT>({0, 4});

    SteinerTreeT algo(G, terminals);
    algo.run();

    EXPECT_EQ(algo.getTotalWeight(), typename SteinerTreeT::WeightSumT{2});
    EXPECT_FALSE(algo.getTree().hasNode(1));
    EXPECT_FALSE(algo.getTree().hasNode(3));
    TestFixture::expectValidTree(G, algo.getTree(), terminals);
}

TYPED_TEST_P(SteinerTreeGTest, testRunTwice) {
    using SteinerTreeT = typename TestFixture::SteinerTreeT;

    typename TestFixture::GraphT G(3);
    G.addEdge(0, 1);
    G.addEdge(1, 2);

    SteinerTreeT algo(G, std::to_array<typename TestFixture::NodeT>({0, 2}));
    algo.run();
    algo.run();

    EXPECT_EQ(algo.getTotalWeight(), typename SteinerTreeT::WeightSumT{2});
    EXPECT_EQ(algo.getTree().numberOfEdges(), 2);
}

TYPED_TEST_P(SteinerTreeGTest, testInvalidInput) {
    using GraphT = typename TestFixture::GraphT;
    using NodeT = typename TestFixture::NodeT;
    using SteinerTreeT = typename TestFixture::SteinerTreeT;

    EXPECT_THROW(SteinerTreeT(GraphT(2, false, true), std::to_array<NodeT>({0})),
                 std::runtime_error);
    EXPECT_THROW(SteinerTreeT(GraphT(2), std::to_array<NodeT>({5})), std::runtime_error);

    if constexpr (std::is_signed_v<typename TestFixture::EdgeWeightT>) {
        GraphT negative(2, true);
        negative.addEdge(0, 1, -1);
        EXPECT_THROW(SteinerTreeT(negative, std::to_array<NodeT>({0})), std::runtime_error);
    }
}

TYPED_TEST_P(SteinerTreeGTest, testResultsBeforeRunThrow) {

    typename TestFixture::GraphT G(2);
    typename TestFixture::SteinerTreeT algo(G, std::to_array<typename TestFixture::NodeT>({0}));
    EXPECT_THROW(algo.getTree(), std::runtime_error);
    EXPECT_THROW(algo.getTotalWeight(), std::runtime_error);
    EXPECT_THROW(algo.hasTree(), std::runtime_error);
}

REGISTER_TYPED_TEST_SUITE_P(SteinerTreeGTest, testUsesASteinerNode,
                            testPathBetweenTwoTerminalsIsShortest,
                            testInducedTreeIsImprovedAndPruned, testLongPathsDoNotOverflow,
                            testUnweightedGraph, testAllNodesAreTerminals,
                            testDisconnectedTerminals, testComponentsWithoutTerminalsAreIgnored,
                            testFewTerminals, testZeroWeightsAndSelfLoops, testDeletedNodes,
                            testRunTwice, testInvalidInput, testResultsBeforeRunThrow);

using SteinerTreeTestTypes =
    ::testing::Types<SteinerTreeConfig<node, edgeweight>, SteinerTreeConfig<uint32_t, float>,
                     SteinerTreeConfig<uint16_t, double>, SteinerTreeConfig<uint64_t, uint64_t>,
                     SteinerTreeConfig<int, int>, SteinerTreeConfig<uint32_t, uint8_t>>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestSteinerTree, SteinerTreeGTest, SteinerTreeTestTypes, );

} // namespace
} // namespace NetworKit
