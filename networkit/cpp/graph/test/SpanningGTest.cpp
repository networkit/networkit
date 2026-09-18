/*
 * SpanningGTest.cpp
 *
 *  Created on: 03.09.2015
 *      Author: Henning
 */

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/generators/BarabasiAlbertGenerator.hpp>
#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/graph/KruskalMSF.hpp>
#include <networkit/graph/PrimMSF.hpp>
#include <networkit/graph/RandomMaximumSpanningForest.hpp>
#include <networkit/graph/SpanningForest.hpp>
#include <networkit/graph/UnionMaximumSpanningForest.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class SpanningGTest : public testing::Test {};

// check that each node has an edge in the spanning tree (if it had one before)
template <typename GraphT>
inline void isValidForest(const GraphT &g, const GraphT &t) {
    using NodeT = typename GraphT::NodeT;

    std::vector<NodeT> graphNodes;
    graphNodes.reserve(g.numberOfNodes());
    g.forNodes([&](NodeT u) { graphNodes.push_back(u); });

    std::vector<NodeT> forestNodes;
    forestNodes.reserve(t.numberOfNodes());
    t.forNodes([&](NodeT u) { forestNodes.push_back(u); });

    EXPECT_THAT(forestNodes, testing::UnorderedElementsAreArray(graphNodes));
    t.forNodes([&](NodeT u) { EXPECT_TRUE(t.degree(u) > 0 || g.degree(u) == 0); });
}

template <class NodeT_, class EdgeWeightT_>
struct PrimMSFConfig {
    using NodeT = NodeT_;
    using EdgeWeightT = EdgeWeightT_;
};

template <class TestT>
class PrimMSFGTest : public testing::Test {
public:
    using NodeT = typename TestT::NodeT;
    using EdgeWeightT = typename TestT::EdgeWeightT;
    using GraphT = AdjListGraph<NodeT, EdgeWeightT>;
    using PrimMSFT = GenericPrimMSF<GraphT>;

    GraphT weightedMSTWithUnitWeights() const {
        GraphT g(5, true);
        g.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{1});
        g.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{1});
        g.addEdge(NodeT{1}, NodeT{3}, EdgeWeightT{1});
        g.addEdge(NodeT{3}, NodeT{4}, EdgeWeightT{1});
        g.addEdge(NodeT{1}, NodeT{4}, EdgeWeightT{1});
        return g;
    }

    GraphT weightedMSFWithUnitWeights() const {
        GraphT g(6, true);
        g.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{1});
        g.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{1});
        g.addEdge(NodeT{2}, NodeT{0}, EdgeWeightT{1});
        g.addEdge(NodeT{3}, NodeT{4}, EdgeWeightT{1});
        g.addEdge(NodeT{4}, NodeT{5}, EdgeWeightT{1});
        g.addEdge(NodeT{5}, NodeT{3}, EdgeWeightT{1});
        return g;
    }

    GraphT weightedMSTWithNonUnitWeights() const {
        GraphT g(4, true);
        g.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{1});
        g.addEdge(NodeT{0}, NodeT{2}, EdgeWeightT{1});
        g.addEdge(NodeT{0}, NodeT{3}, EdgeWeightT{1});
        g.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{2});
        g.addEdge(NodeT{2}, NodeT{3}, EdgeWeightT{2});
        return g;
    }

    GraphT weightedMSFWithNonUnitWeights() const {
        GraphT g(6, true);
        g.addEdge(NodeT{0}, NodeT{1}, EdgeWeightT{1});
        g.addEdge(NodeT{1}, NodeT{2}, EdgeWeightT{2});
        g.addEdge(NodeT{2}, NodeT{0}, EdgeWeightT{3});
        g.addEdge(NodeT{3}, NodeT{4}, EdgeWeightT{1});
        g.addEdge(NodeT{4}, NodeT{5}, EdgeWeightT{2});
        g.addEdge(NodeT{5}, NodeT{3}, EdgeWeightT{3});
        return g;
    }

    GraphT unweightedMSF() const {
        GraphT g(6);
        g.addEdge(NodeT{0}, NodeT{1});
        g.addEdge(NodeT{1}, NodeT{2});
        g.addEdge(NodeT{2}, NodeT{0});
        g.addEdge(NodeT{3}, NodeT{4});
        g.addEdge(NodeT{4}, NodeT{5});
        g.addEdge(NodeT{5}, NodeT{3});
        return g;
    }
};

TYPED_TEST_SUITE_P(PrimMSFGTest);

TEST_F(SpanningGTest, testSpanningForest) {
    METISGraphReader reader;
    std::vector<std::string> graphs = {"karate", "jazz", "celegans_metabolic"};

    for (const auto &graphname : graphs) {
        std::string filename = "input/" + graphname + ".graph";
        Graph G = reader.read(filename);
        SpanningForest msf(G);
        msf.run();
        Graph T = msf.getForest();

        INFO("tree / graph edges: ", T.numberOfEdges(), " / ", G.numberOfEdges());

        isValidForest(G, T);
    }
}

TEST_F(SpanningGTest, testRandomMaximumSpanningForest) {
    METISGraphReader reader;
    std::vector<std::string> graphs = {"karate", "jazz", "celegans_metabolic"};

    for (const auto &graphname : graphs) {
        std::string filename = "input/" + graphname + ".graph";
        Graph G = reader.read(filename);

        RandomMaximumSpanningForest rmsf(G);
        rmsf.run();
        Graph T = rmsf.getMSF();

        INFO("tree / graph edges: ", T.numberOfEdges(), " / ", G.numberOfEdges());

        // check that each node has an edge in the spanning tree (if it had one before)
        T.forNodes([&](node u) { EXPECT_TRUE(T.degree(u) > 0 || G.degree(u) == 0); });
        T.forEdges([&](node u, node v) { EXPECT_TRUE(rmsf.inMSF(u, v)); });
    }
}

TEST_F(SpanningGTest, testUnionMaximumSpanningForest) {
    METISGraphReader reader;
    std::vector<std::string> graphs = {"karate", "jazz", "celegans_metabolic"};

    for (const auto &graphname : graphs) {
        std::string filename = "input/" + graphname + ".graph";
        Graph G = reader.read(filename);

        UnionMaximumSpanningForest umsf(G);
        umsf.run();
        Graph T = umsf.getUMSF();

        INFO("tree / graph edges: ", T.numberOfEdges(), " / ", G.numberOfEdges());

        // check that each node has an edge in the spanning tree (if it had one before)
        T.forNodes([&](node u) { EXPECT_TRUE(T.degree(u) > 0 || G.degree(u) == 0); });
        T.forEdges([&](node u, node v) { EXPECT_TRUE(umsf.inUMSF(u, v)); });
    }
}

TEST_F(SpanningGTest, testKruskalMinSpanningForest) {
    METISGraphReader reader;
    std::vector<std::string> graphs = {"karate", "jazz", "celegans_metabolic"};

    for (const auto &graphname : graphs) {
        std::string filename = "input/" + graphname + ".graph";
        Graph G = reader.read(filename);
        KruskalMSF msf(G);
        msf.run();
        Graph T = msf.getForest();

        isValidForest(G, T);
    }
}

TEST_F(SpanningGTest, testKruskalMinimumSpanningForestUnweightedGraph) {
    Graph g(5, false);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(1, 3);
    g.addEdge(3, 4);
    g.addEdge(1, 4);
    g.indexEdges();

    KruskalMSF msf(g);
    msf.run();

    isValidForest(g, msf.getForest());
    EXPECT_EQ(msf.getTotalWeight(), 4);
}

TEST_F(SpanningGTest, testKruskalMinimumSpanningForestIsMSTUnitWeights) {
    Graph g(5, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 1);
    g.addEdge(3, 4, 1);
    g.addEdge(1, 4, 1);
    g.indexEdges();

    KruskalMSF msf(g);
    msf.run();
    Graph T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 4);
}

TEST_F(SpanningGTest, testKruskalMinimumSpanningForestIsMSFUnitWeights) {
    Graph g(6, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(2, 0, 1);
    g.addEdge(3, 4, 1);
    g.addEdge(4, 5, 1);
    g.addEdge(5, 3, 1);
    g.indexEdges();

    KruskalMSF msf(g);
    msf.run();
    Graph T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 4);
}

TEST_F(SpanningGTest, testKruskalMinimumSpanningForestIsMSTNonUnitWeights) {
    Graph g(4, true);
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 1);
    g.addEdge(0, 3, 1);
    g.addEdge(1, 2, 2);
    g.addEdge(2, 3, 2);
    g.indexEdges();

    KruskalMSF msf(g);
    msf.run();
    Graph T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 3);
}

TEST_F(SpanningGTest, testKruskalMinimumSpanningForestIsMSFNonUnitWeights) {
    Graph g(6, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 2);
    g.addEdge(2, 0, 3);
    g.addEdge(3, 4, 1);
    g.addEdge(4, 5, 2);
    g.addEdge(5, 3, 3);
    g.indexEdges();

    KruskalMSF msf(g);
    msf.run();
    Graph T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 6);
}

TYPED_TEST_P(PrimMSFGTest, testThrowsForDirectedGraph) {
    typename TestFixture::GraphT g(5, true, true);
    try {
        typename TestFixture::PrimMSFT msf(g);
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "The graph is not an undirected graph.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(SpanningGTest, testPrimMinSpanningForest) {
    METISGraphReader reader;
    std::vector<std::string> graphs = {"karate", "jazz", "celegans_metabolic"};

    for (const auto &graphname : graphs) {
        std::string filename = "input/" + graphname + ".graph";
        Graph G = reader.read(filename);
        PrimMSF msf(G);
        msf.run();
        Graph T = msf.getForest();

        isValidForest(G, T);
    }
}

TYPED_TEST_P(PrimMSFGTest, testMinimumSpanningForestIsMSTUnitWeights) {
    auto g = this->weightedMSTWithUnitWeights();

    typename TestFixture::PrimMSFT msf(g);
    msf.run();
    const auto &T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 4);
}

TYPED_TEST_P(PrimMSFGTest, testMinimumSpanningForestIsMSFUnitWeights) {
    auto g = this->weightedMSFWithUnitWeights();

    typename TestFixture::PrimMSFT msf(g);
    msf.run();
    const auto &T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 4);
}

TYPED_TEST_P(PrimMSFGTest, testMinimumSpanningForestIsMSTNonUnitWeights) {
    auto g = this->weightedMSTWithNonUnitWeights();

    typename TestFixture::PrimMSFT msf(g);
    msf.run();
    const auto &T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 3);
}

TYPED_TEST_P(PrimMSFGTest, testMinimumSpanningForestIsMSFNonUnitWeights) {
    auto g = this->weightedMSFWithNonUnitWeights();

    typename TestFixture::PrimMSFT msf(g);
    msf.run();
    const auto &T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 6);
}

TYPED_TEST_P(PrimMSFGTest, testMinimumSpanningForestIsMSFUnweighted) {
    auto g = this->unweightedMSF();

    typename TestFixture::PrimMSFT msf(g);
    msf.run();
    const auto &T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 4);
}

REGISTER_TYPED_TEST_SUITE_P(PrimMSFGTest, testThrowsForDirectedGraph,
                            testMinimumSpanningForestIsMSTUnitWeights,
                            testMinimumSpanningForestIsMSFUnitWeights,
                            testMinimumSpanningForestIsMSTNonUnitWeights,
                            testMinimumSpanningForestIsMSFNonUnitWeights,
                            testMinimumSpanningForestIsMSFUnweighted);

using PrimMSFTestTypes = ::testing::Types<PrimMSFConfig<node, edgeweight>,
                                          PrimMSFConfig<uint32_t, float>, PrimMSFConfig<int, int>>;

INSTANTIATE_TYPED_TEST_SUITE_P(TestPrimMSF, PrimMSFGTest, PrimMSFTestTypes, );

} /* namespace NetworKit */
