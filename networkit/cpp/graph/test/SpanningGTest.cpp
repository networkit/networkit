/*
 * SpanningGTest.cpp
 *
 *  Created on: 03.09.2015
 *      Author: Henning
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/graph/KruskalMSF.hpp>
#include <networkit/graph/PrimMSF.hpp>
#include <networkit/graph/RandomMaximumSpanningForest.hpp>
#include <networkit/graph/SpanningForest.hpp>
#include <networkit/graph/UnionMaximumSpanningForest.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class SpanningGTest : public testing::Test {};

// check that each node has an edge in the spanning tree (if it had one before)
inline void isValidForest(const Graph &g, const Graph &t) {
    t.forNodes([&](node u) { EXPECT_TRUE(t.degree(u) > 0 || g.degree(u) == 0); });
}

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

TEST_F(SpanningGTest, testKruskalMinimumSpanningForestIsMSTUnitWeights) {
    GraphW g(5, true);
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
    GraphW g(6, true);
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
    GraphW g(4, true);
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
    GraphW g(6, true);
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

TEST_F(SpanningGTest, testPrimThrowsForDirectedGraph) {
    GraphW g(5, true, true);
    try {
        PrimMSF msf(g);
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

TEST_F(SpanningGTest, testPrimMinimumSpanningForestIsMSTUnitWeights) {
    GraphW g(5, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 1);
    g.addEdge(3, 4, 1);
    g.addEdge(1, 4, 1);
    g.indexEdges();

    PrimMSF msf(g);
    msf.run();
    Graph T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 4);
}

TEST_F(SpanningGTest, testPrimMinimumSpanningForestIsMSFUnitWeights) {
    GraphW g(6, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 1);
    g.addEdge(2, 0, 1);
    g.addEdge(3, 4, 1);
    g.addEdge(4, 5, 1);
    g.addEdge(5, 3, 1);
    g.indexEdges();

    PrimMSF msf(g);
    msf.run();
    Graph T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 4);
}

TEST_F(SpanningGTest, testPrimMinimumSpanningForestIsMSTNonUnitWeights) {
    GraphW g(4, true);
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 1);
    g.addEdge(0, 3, 1);
    g.addEdge(1, 2, 2);
    g.addEdge(2, 3, 2);
    g.indexEdges();

    PrimMSF msf(g);
    msf.run();
    Graph T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 3);
}

TEST_F(SpanningGTest, testPrimMinimumSpanningForestIsMSFNonUnitWeights) {
    GraphW g(6, true);
    g.addEdge(0, 1, 1);
    g.addEdge(1, 2, 2);
    g.addEdge(2, 0, 3);
    g.addEdge(3, 4, 1);
    g.addEdge(4, 5, 2);
    g.addEdge(5, 3, 3);
    g.indexEdges();

    PrimMSF msf(g);
    msf.run();
    Graph T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 6);
}

TEST_F(SpanningGTest, testPrimMinimumSpanningForestIsMSFUnweighted) {
    GraphW g(6);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 0);
    g.addEdge(3, 4);
    g.addEdge(4, 5);
    g.addEdge(5, 3);
    g.indexEdges();

    PrimMSF msf(g);
    msf.run();
    Graph T = msf.getForest();

    isValidForest(g, T);
    EXPECT_EQ(msf.getTotalWeight(), 4);
}

} /* namespace NetworKit */
