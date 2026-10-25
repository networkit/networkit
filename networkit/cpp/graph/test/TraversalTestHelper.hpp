/*  TraversalTestHelper.cpp
 *
 *  Created on: 25.10.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_GRAPH_TEST_TRAVERSAL_TEST_HELPER_HPP_
#define NETWORKIT_GRAPH_TEST_TRAVERSAL_TEST_HELPER_HPP_

#include <cstdint>

#include <gtest/gtest.h>

#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

template <class GraphT_>
struct TraversalConfig {
    using GraphT = GraphT_;
};

template <class TestT>
class GenericTraversalGTest : public testing::Test {
public:
    using GraphT = typename TestT::GraphT;
    using NodeT = typename GraphT::NodeT;
    using EdgeWeightT = typename GraphT::EdgeWeightT;

    GraphT unweightedGraph() const {
        GraphT G(7, false);
        G.addEdge(NodeT{0}, NodeT{1});
        G.addEdge(NodeT{0}, NodeT{2});
        G.addEdge(NodeT{1}, NodeT{3});
        G.addEdge(NodeT{2}, NodeT{4});
        G.addEdge(NodeT{3}, NodeT{5});
        G.addEdge(NodeT{4}, NodeT{6});
        return G;
    }

    GraphT weightedGraph() const {
        GraphT G(7, true);
        G.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{10});
        G.addEdge(NodeT{0}, NodeT{2}, EdgeWeightT{20});
        G.addEdge(NodeT{1}, NodeT{3}, EdgeWeightT{30});
        G.addEdge(NodeT{2}, NodeT{4}, EdgeWeightT{40});
        G.addEdge(NodeT{3}, NodeT{5}, EdgeWeightT{50});
        G.addEdge(NodeT{4}, NodeT{6}, EdgeWeightT{60});
        G.addEdge(NodeT{2}, NodeT{3}, EdgeWeightT{70});
        return G;
    }

    GraphT pathGraph() const {
        GraphT G(5, false);
        G.addEdge(NodeT{0}, NodeT{1});
        G.addEdge(NodeT{1}, NodeT{2});
        G.addEdge(NodeT{2}, NodeT{3});
        G.addEdge(NodeT{3}, NodeT{4});
        return G;
    }

    GraphT directedGraph() const {
        GraphT G(6, false, true);
        G.addEdge(NodeT{0}, NodeT{1});
        G.addEdge(NodeT{1}, NodeT{2});
        G.addEdge(NodeT{2}, NodeT{4});
        G.addEdge(NodeT{3}, NodeT{0});
        G.addEdge(NodeT{4}, NodeT{4});
        return G;
    }
};

using TraversalTestTypes = ::testing::Types<
    TraversalConfig<Graph>, TraversalConfig<AdjListGraph<uint32_t, float>>,
    TraversalConfig<AdjListGraph<uint16_t, double>>,
    TraversalConfig<AdjListGraph<uint64_t, uint64_t>>, TraversalConfig<AdjListGraph<int, int>>,
    TraversalConfig<AdjListGraph<int32_t, double>>, TraversalConfig<AdjListGraph<int64_t, float>>>;

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_TEST_TRAVERSAL_TEST_HELPER_HPP_
