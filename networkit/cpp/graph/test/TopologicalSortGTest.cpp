/*
 * TopologicalSort.hpp
 *
 *  Created on: 22.11.2021
 *      Author: Fabian Brandt-Tumescheit
 */

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/graph/TopologicalSort.hpp>

namespace NetworKit {
namespace {

template <class NodeT_, class EdgeWeightT_>
struct TopologicalSortConfig {
    using NodeT = NodeT_;
    using EdgeWeightT = EdgeWeightT_;
};

template <class TestT>
class TopologicalSortGTest : public testing::Test {
public:
    using NodeT = typename TestT::NodeT;
    using EdgeWeightT = typename TestT::EdgeWeightT;
    using GraphT = AdjListGraph<NodeT, EdgeWeightT>;
    using TopologicalSortT = GenericTopologicalSort<GraphT>;
    using NodeIdMapping = typename TopologicalSortT::NodeIdMapping;

    GraphT inputGraph(bool directed) const {
        GraphT G(5, false, directed);

        /**
         * /--> 1 --> 3
         * |    ^
         * 0    |
         * |    |
         * \--> 2 <-- 4
         */

        G.addEdge(NodeT{0}, NodeT{1});
        G.addEdge(NodeT{0}, NodeT{2});
        G.addEdge(NodeT{2}, NodeT{1});
        G.addEdge(NodeT{1}, NodeT{3});
        G.addEdge(NodeT{4}, NodeT{2});

        return G;
    }

    NodeIdMapping makeMapping() const {
        NodeIdMapping mapping;
        mapping[NodeT{0}] = 0;
        mapping[NodeT{1}] = 1;
        mapping[NodeT{2}] = 2;
        mapping[NodeT{4}] = 3;

        return mapping;
    }

    void assertTopological(const GraphT &G, const std::vector<NodeT> &sort) const {
        std::vector<NodeT> nodes;
        nodes.reserve(G.numberOfNodes());
        G.forNodes([&](NodeT u) { nodes.push_back(u); });

        std::unordered_map<NodeT, index> indices;
        EXPECT_EQ(sort.size(), G.numberOfNodes());
        EXPECT_THAT(sort, testing::UnorderedElementsAreArray(nodes));
        G.forNodes([&](NodeT u) {
            const auto it = std::find(sort.begin(), sort.end(), u);
            EXPECT_NE(it, sort.end());
            indices[u] = std::distance(sort.begin(), it);
        });
        G.forNodes([&](NodeT u) {
            G.forNeighborsOf(u, [&](NodeT v) { EXPECT_LT(indices[u], indices[v]); });
        });
    }
};

TYPED_TEST_SUITE_P(TopologicalSortGTest);

TYPED_TEST_P(TopologicalSortGTest, testTopologicalSort) {
    auto G = this->inputGraph(true);

    typename TestFixture::TopologicalSortT topSort(G);
    topSort.run();
    const auto &res = topSort.getResult();

    this->assertTopological(G, res);
}

TYPED_TEST_P(TopologicalSortGTest, testRepeatedRuns) {
    auto G = this->inputGraph(true);

    typename TestFixture::TopologicalSortT topSort(G);
    topSort.run();
    const auto res = topSort.getResult();
    topSort.run();
    const auto res2 = topSort.getResult();

    EXPECT_THAT(res2, testing::ElementsAreArray(res));
}

TYPED_TEST_P(TopologicalSortGTest, testRejectGraphWithCycles) {
    using NodeT = typename TestFixture::NodeT;

    auto G = this->inputGraph(true);
    G.addEdge(NodeT{3}, NodeT{4});

    typename TestFixture::TopologicalSortT topSort(G);
    EXPECT_THROW(topSort.run(), std::runtime_error);
}

TYPED_TEST_P(TopologicalSortGTest, testRejectUndirectedGraph) {
    EXPECT_THROW(typename TestFixture::TopologicalSortT(this->inputGraph(false)),
                 std::runtime_error);
}

TYPED_TEST_P(TopologicalSortGTest, testNonContinuousNodeIds) {
    using NodeT = typename TestFixture::NodeT;

    auto G = this->inputGraph(true);
    G.removeNode(NodeT{3});

    typename TestFixture::TopologicalSortT topSort(G);
    topSort.run();
    const auto &res = topSort.getResult();

    this->assertTopological(G, res);
}

TYPED_TEST_P(TopologicalSortGTest, testCustomNodeIdMapping) {
    using NodeT = typename TestFixture::NodeT;

    auto G = this->inputGraph(true);
    G.removeNode(NodeT{3});
    auto mapping = this->makeMapping();

    typename TestFixture::TopologicalSortT topSort(G, mapping);
    topSort.run();
    const auto &res = topSort.getResult();

    this->assertTopological(G, res);
}

TYPED_TEST_P(TopologicalSortGTest, testWrongSizeOfNodeIdMapping) {
    using NodeT = typename TestFixture::NodeT;

    auto G = this->inputGraph(true);
    G.removeNode(NodeT{3});
    auto mapping = this->makeMapping();
    mapping[NodeT{5}] = 4;

    EXPECT_THROW(typename TestFixture::TopologicalSortT(G, mapping), std::runtime_error);
}

TYPED_TEST_P(TopologicalSortGTest, testNonContinuousNodeIdMapping) {
    using NodeT = typename TestFixture::NodeT;

    auto G = this->inputGraph(true);
    G.removeNode(NodeT{3});
    auto mapping = this->makeMapping();
    mapping[NodeT{1}] = 4;
    EXPECT_THROW(typename TestFixture::TopologicalSortT(G, mapping, true), std::runtime_error);

    mapping = this->makeMapping();
    mapping.erase(NodeT{1});
    // to get correct size
    mapping[NodeT{5}] = 5;
    EXPECT_THROW(typename TestFixture::TopologicalSortT(G, mapping, true), std::runtime_error);

    typename TestFixture::TopologicalSortT topSort(G, mapping);
    EXPECT_THROW(topSort.run(), std::runtime_error);
}

TYPED_TEST_P(TopologicalSortGTest, testNonInjectiveNodeIdMapping) {
    using NodeT = typename TestFixture::NodeT;

    auto G = this->inputGraph(true);
    G.removeNode(NodeT{3});
    auto mapping = this->makeMapping();
    mapping[NodeT{2}] = 1;

    EXPECT_THROW(typename TestFixture::TopologicalSortT(G, mapping, true), std::runtime_error);
}

REGISTER_TYPED_TEST_SUITE_P(TopologicalSortGTest, testTopologicalSort, testRepeatedRuns,
                            testRejectGraphWithCycles, testRejectUndirectedGraph,
                            testNonContinuousNodeIds, testCustomNodeIdMapping,
                            testWrongSizeOfNodeIdMapping, testNonContinuousNodeIdMapping,
                            testNonInjectiveNodeIdMapping);

using TopologicalSortTestTypes =
    ::testing::Types<TopologicalSortConfig<node, edgeweight>,
                     TopologicalSortConfig<uint32_t, float>, TopologicalSortConfig<int, int>>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestTopologicalSort, TopologicalSortGTest,
                               TopologicalSortTestTypes, );

} // namespace
} // namespace NetworKit
