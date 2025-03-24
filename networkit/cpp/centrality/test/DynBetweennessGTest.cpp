/*
 * DynBetweennessGTest.cpp
 *
 *  Created on: 05.08.2014
 *      Author: ebergamini, cls
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/centrality/ApproxBetweenness.hpp>
#include <networkit/centrality/Betweenness.hpp>
#include <networkit/centrality/DynApproxBetweenness.hpp>
#include <networkit/centrality/DynBetweenness.hpp>
#include <networkit/centrality/DynBetweennessOneNode.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class DynBetweennessGTest : public testing::TestWithParam<std::pair<bool, bool>> {
protected:
    bool isDirected() const noexcept { return GetParam().first; };
    bool isWeighted() const noexcept { return GetParam().second; };
    Graph generateSmallGraph() const;

    static constexpr double epsilon = 0.1, delta = 0.1;

    void compareAgainstBaseline(const Graph &G, const std::vector<double> &apxScores,
                                const std::vector<double> &exactScores, double normalized = false,
                                double err = epsilon) const {
        const auto n = static_cast<double>(G.numberOfNodes());
        const auto normFactor = normalized ? n * (n - 1) : 1.;
        G.forNodes([&](node u) { EXPECT_NEAR(apxScores[u], exactScores[u] / normFactor, err); });
    }

    std::pair<node, node> getNonAdjacentNodes(const Graph &G) const {
        node u, v;
        do {
            u = GraphTools::randomNode(G), v = GraphTools::randomNode(G);
        } while (u == v || G.hasEdge(u, v));
        return {u, v};
    }
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, DynBetweennessGTest,
                         testing::Values(std::make_pair(false, false), std::make_pair(true, false),
                                         std::make_pair(false, true), std::make_pair(true, true)));

Graph DynBetweennessGTest::generateSmallGraph() const {
    /* Graph:
       0    3   6
        \  / \ /
         2    5
        /  \ / \
       1    4   7
    */
    Graph G(8, isWeighted(), isDirected());

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.addEdge(5, 6);
    G.addEdge(5, 7);

    if (isDirected()) {
        G.addEdge(4, 1);
        G.addEdge(3, 0);
        G.addEdge(5, 2);
    }

    if (isWeighted()) {
        Aux::Random::setSeed(42, false);
        GraphTools::randomizeWeights(G);
    }
    return G;
}

TEST_P(DynBetweennessGTest, runDynApproxBetweennessSmallGraph) {
    Graph G = generateSmallGraph();
    DynApproxBetweenness dynbc(G, epsilon, delta);
    Betweenness bc(G);
    dynbc.run();
    bc.run();
    compareAgainstBaseline(G, dynbc.scores(), bc.scores(), true);

    GraphEvent event(GraphEvent::EDGE_ADDITION, 0, 6);
    G.addEdge(event.u, event.v);
    bc.run();
    dynbc.update(event);
    compareAgainstBaseline(G, dynbc.scores(), bc.scores(), true);
}

TEST_P(DynBetweennessGTest, runDynApproxBetweennessSmallGraphEdgeDeletion) {
    Graph G = generateSmallGraph();
    DynApproxBetweenness dynbc(G, epsilon, delta);
    Betweenness bc(G);
    dynbc.run();
    bc.run();
    compareAgainstBaseline(G, dynbc.scores(), bc.scores(), true);

    GraphEvent event(GraphEvent::EDGE_REMOVAL, 3, 5);
    G.removeEdge(event.u, event.v);
    bc.run();
    dynbc.update(event);
    compareAgainstBaseline(G, dynbc.scores(), bc.scores(), true);
}

TEST_P(DynBetweennessGTest, testDynApproxBetweenessGeneratedGraph) {
    Aux::Random::setSeed(42, false);
    ErdosRenyiGenerator generator(100, 0.25, isDirected());
    Graph G = generator.generate();
    if (isWeighted()) {
        GraphTools::randomizeWeights(G);
    }

    DynApproxBetweenness dynbc(G, epsilon, delta);
    Betweenness bc(G);
    dynbc.run();
    for (int i = 0; i < 10; ++i) {
        auto [v1, v2] = getNonAdjacentNodes(G);
        G.addEdge(v1, v2);

        GraphEvent event(GraphEvent::EDGE_ADDITION, v1, v2);
        dynbc.update(event);
        bc.run();
        compareAgainstBaseline(G, dynbc.scores(), bc.scores(), true);
    }
}

TEST_P(DynBetweennessGTest, runDynApproxBetweenessGeneratedGraphEdgeDeletion) {
    Aux::Random::setSeed(42, false);
    ErdosRenyiGenerator generator(100, 0.25, isDirected());
    Graph G = generator.generate();
    if (isWeighted()) {
        Aux::Random::setSeed(42, false);
        GraphTools::randomizeWeights(G);
    }

    DynApproxBetweenness dynbc(G, epsilon, delta);
    Betweenness bc(G);
    dynbc.run();
    for (count i = 0; i < 10; i++) {
        auto randomEdge = GraphTools::randomEdge(G);
        G.removeEdge(randomEdge.first, randomEdge.second);
        GraphEvent event(GraphEvent::EDGE_REMOVAL, randomEdge.first, randomEdge.second);
        dynbc.update(event);
        bc.run();
        compareAgainstBaseline(G, dynbc.scores(), bc.scores(), true);
    }
}

TEST_F(DynBetweennessGTest, runDynVsStatic) {
    Graph G = METISGraphReader{}.read("input/celegans_metabolic.graph");

    DynApproxBetweenness dynbc(G, epsilon, delta, false);
    ApproxBetweenness bc(G, epsilon, delta);
    dynbc.run();
    bc.run();
    compareAgainstBaseline(G, dynbc.scores(), bc.scores());

    std::vector<GraphEvent> batch;
    for (int i = 0; i < 10; ++i) {
        auto [u, v] = getNonAdjacentNodes(G);
        G.addEdge(u, v);
        batch.emplace_back(GraphEvent::EDGE_ADDITION, u, v);
    }
    bc.run();
    dynbc.updateBatch(batch);
    compareAgainstBaseline(G, dynbc.scores(), bc.scores());
}

TEST_F(DynBetweennessGTest, runDynVsStaticEdgeDeletion) {
    Graph G = METISGraphReader{}.read("input/celegans_metabolic.graph");

    DynApproxBetweenness dynbc(G, epsilon, delta, false);
    ApproxBetweenness bc(G, epsilon, delta);
    dynbc.run();
    bc.run();
    compareAgainstBaseline(G, dynbc.scores(), bc.scores());

    std::vector<GraphEvent> batch;
    for (int i = 0; i < 10; ++i) {
        auto randomEdge = GraphTools::randomEdge(G);
        G.removeEdge(randomEdge.first, randomEdge.second);
        batch.emplace_back(GraphEvent::EDGE_REMOVAL, randomEdge.first, randomEdge.second);
    }
    bc.run();
    dynbc.updateBatch(batch);
    compareAgainstBaseline(G, dynbc.scores(), bc.scores());
}

TEST_F(DynBetweennessGTest, runDynVsStaticCaseInsertDirected) {
    Aux::Random::setSeed(0, false);

    for (count n = 2; n <= 25; n++)
        for (count t = 0; t < 100; t++) {
            auto g = ErdosRenyiGenerator(n, 0.3, true).generate();

            if (g.numberOfEdges() == g.numberOfNodes() * (g.numberOfNodes() - 1))
                continue;

            auto [u, v] = getNonAdjacentNodes(g);
            auto ibet = DynBetweenness(g);

            ibet.run();

            auto ne = GraphEvent(GraphEvent::EDGE_ADDITION, u, v);
            g.addEdge(u, v);
            ibet.update(ne);

            EXPECT_TRUE(g.hasEdge(u, v));
            auto brandes = Betweenness(g);
            brandes.run();
            g.forNodes([&](node w) { EXPECT_NEAR(brandes.score(w), ibet.score(w), 1e-8); });
        }
}

TEST_F(DynBetweennessGTest, runDynVsStaticCaseInsertUndirected) {
    Aux::Random::setSeed(0, false);

    for (count n = 2; n <= 25; n++)
        for (count t = 0; t < 100; t++) {
            auto g = ErdosRenyiGenerator(n, 0.3, false).generate();
            if (g.numberOfEdges() == g.numberOfNodes() * (g.numberOfNodes() - 1) / 2)
                continue;

            auto [u, v] = getNonAdjacentNodes(g);
            auto ibet = DynBetweenness(g);

            ibet.run();

            auto ne = GraphEvent(GraphEvent::EDGE_ADDITION, u, v);
            g.addEdge(u, v);
            ibet.update(ne);
            EXPECT_TRUE(g.hasEdge(u, v));

            auto brandes = Betweenness(g);
            brandes.run();
            g.forNodes([&](node w) { EXPECT_NEAR(brandes.score(w), ibet.score(w), 1e-8); });
        }
}

TEST_P(DynBetweennessGTest, testDynamicBetweennessOneNode) {
    // for each of the 8 focus nodes in the small graph
    for (node x = 0; x < 8; ++x) {
        Graph G = generateSmallGraph();

        DynBetweennessOneNode dynb(G, x);
        dynb.run();

        Betweenness b(G);
        b.run();

        EXPECT_NEAR(dynb.getbcx(), b.score(x), 0.01);

        for (count i = 0; i < 10; ++i) {
            auto e = getNonAdjacentNodes(G);
            edgeweight ew = 1.0;
            if (G.isWeighted()) {
                ew = Aux::Random::real(1);
            }
            G.addEdge(e.first, e.second, ew);
            dynb.update(GraphEvent(GraphEvent::EDGE_ADDITION, e.first, e.second, ew));
        }

        Betweenness b_updated(G);
        b_updated.run();

        EXPECT_NEAR(dynb.getbcx(), b_updated.score(x), 0.01);
    }
}

} /* namespace NetworKit */
