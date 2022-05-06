// no-networkit-format
/*
 * DynBetweennessGTest.cpp
 *
 *  Created on: 05.08.2014
 *      Author: ebergamini, cls
 */

#include <gtest/gtest.h>

#include <networkit/centrality/Betweenness.hpp>
#include <networkit/centrality/DynApproxBetweenness.hpp>
#include <networkit/centrality/ApproxBetweenness.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/generators/DorogovtsevMendesGenerator.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/centrality/DynBetweenness.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

class DynBetweennessGTest: public testing::TestWithParam<std::pair<bool, bool>> {
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
                        std::make_pair(false, true),
                        std::make_pair(true, true)));

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
        G.forEdges([&](node u, node v) { G.setWeight(u, v, Aux::Random::probability()); });
    }
    return G;
}

TEST_P(DynBetweennessGTest, runDynApproxBetweennessSmallGraph) {
    Graph G = generateSmallGraph();
    DynApproxBetweenness dynbc(G, epsilon, delta);
    Betweenness bc(G);
    dynbc.run();
    bc.run();
    std::vector<double> dynbc_scores = dynbc.scores();
    std::vector<double> bc_scores = bc.scores();
    compareAgainstBaseline(G, dynbc_scores, bc_scores, true);

    GraphEvent event(GraphEvent::EDGE_ADDITION, 0, 6);
    G.addEdge(event.u, event.v);
    bc.run();
    dynbc.update(event);
    dynbc_scores = dynbc.scores();
    bc_scores = bc.scores();
    compareAgainstBaseline(G, dynbc_scores, bc_scores, true);
}

TEST_P(DynBetweennessGTest, runDynApproxBetweennessSmallGraphEdgeDeletion) {
    Graph G = generateSmallGraph();
    DynApproxBetweenness dynbc(G, epsilon, delta);
    Betweenness bc(G);
    dynbc.run();
    bc.run();
    std::vector<double> dynbc_scores = dynbc.scores();
    std::vector<double> bc_scores = bc.scores();
    compareAgainstBaseline(G, dynbc_scores, bc_scores, true);

    GraphEvent event(GraphEvent::EDGE_REMOVAL, 3, 5);
    G.removeEdge(event.u, event.v);
    bc.run();
    dynbc.update(event);
    dynbc_scores = dynbc.scores();
    bc_scores = bc.scores();
    compareAgainstBaseline(G, dynbc_scores, bc_scores, true);
}

TEST_P(DynBetweennessGTest, testDynApproxBetweenessGeneratedGraph) {
    ErdosRenyiGenerator generator(100, 0.25, isDirected());
    Graph G = generator.generate();
    if (isWeighted()) {
        G = GraphTools::toWeighted(G);
        Aux::Random::setSeed(42, false);
        G.forEdges([&](node u, node v) { G.setWeight(u, v, Aux::Random::probability()); });
    }

    DEBUG("Generated graph of dimension ", G.upperNodeIdBound());
    DynApproxBetweenness dynbc(G, epsilon, delta);
    Betweenness bc(G);
    dynbc.run();
    DEBUG("Before the edge insertion: ");
    for (int i = 0; i < 10; ++i) {
        auto [v1, v2] = getNonAdjacentNodes(G);
        G.addEdge(v1, v2);

        GraphEvent event(GraphEvent::EDGE_ADDITION, v1, v2);
        dynbc.update(event);
        bc.run();
        std::vector<double> dynbc_scores = dynbc.scores();
        std::vector<double> bc_scores = bc.scores();
        compareAgainstBaseline(G, dynbc_scores, bc_scores, true);
    }
}

TEST_P(DynBetweennessGTest, runDynApproxBetweenessGeneratedGraphEdgeDeletion) {
    ErdosRenyiGenerator generator(100, 0.25, isDirected());
    Graph G = generator.generate();

    DEBUG("Generated graph of dimension ", G.upperNodeIdBound());
    DynApproxBetweenness dynbc(G, epsilon, delta);
    Betweenness bc(G);
    dynbc.run();
    DEBUG("Before the edge deletion: ");
    for(count i = 0; i < 10; i++) {
        DEBUG("Selecting a random edge");
        auto randomEdge = GraphTools::randomEdge(G);
        DEBUG("Deleting edge number ", i);
        G.removeEdge(randomEdge.first, randomEdge.second);
        GraphEvent event(GraphEvent::EDGE_REMOVAL, randomEdge.first, randomEdge.second);
        dynbc.update(event);
        bc.run();
        std::vector<double> dynbc_scores = dynbc.scores();
        std::vector<double> bc_scores = bc.scores();
        compareAgainstBaseline(G, dynbc_scores, bc_scores, true);
    }
}

TEST_F(DynBetweennessGTest, runDynVsStatic) {
    Graph G = METISGraphReader{}.read("input/celegans_metabolic.graph");

    DEBUG("Initializing DynApproxBetweenness");
    DynApproxBetweenness dynbc(G, epsilon, delta, false);
    DEBUG("Initializing ApproxBetweenness");
    ApproxBetweenness bc(G, epsilon, delta);
    DEBUG("Running DynApproxBetweenness");
    dynbc.run();
    DEBUG("Running ApproxBetweenness");
    bc.run();
    std::vector<double> dynbc_scores = dynbc.scores();
    std::vector<double> bc_scores = bc.scores();
    compareAgainstBaseline(G, dynbc_scores, bc_scores);

    DEBUG("Before the edge insertion: ");
    std::vector<GraphEvent> batch;
    for (int i = 0; i < 10; ++i) {
        auto [u, v] = getNonAdjacentNodes(G);
        G.addEdge(u, v);
        batch.emplace_back(GraphEvent::EDGE_ADDITION, u, v);
    }
    DEBUG("Running ApproxBetweenness (again)");
    bc.run();
    DEBUG("Updating DynApproxBetweenness");
    dynbc.updateBatch(batch);
    DEBUG("Calling DynApproxBetweenness Scores");
    dynbc_scores = dynbc.scores();
    DEBUG("Calling ApproxBetweenness Scores");
    bc_scores = bc.scores();
    compareAgainstBaseline(G, dynbc_scores, bc_scores);
}

TEST_F(DynBetweennessGTest, runDynVsStaticEdgeDeletion) {
    Graph G = METISGraphReader{}.read("input/celegans_metabolic.graph");

    DynApproxBetweenness dynbc(G, epsilon, delta, false);
    ApproxBetweenness bc(G, epsilon, delta);
    dynbc.run();
    bc.run();
    std::vector<double> dynbc_scores = dynbc.scores();
    std::vector<double> bc_scores = bc.scores();
    compareAgainstBaseline(G, dynbc_scores, bc_scores);

    std::vector<GraphEvent> batch;
    for (int i = 0; i < 10; ++i) {
        auto randomEdge = GraphTools::randomEdge(G);
        G.removeEdge(randomEdge.first, randomEdge.second);
        batch.emplace_back(GraphEvent::EDGE_REMOVAL, randomEdge.first, randomEdge.second);
    }
    bc.run();
    dynbc.updateBatch(batch);
    dynbc_scores = dynbc.scores();
    bc_scores = bc.scores();
    compareAgainstBaseline(G, dynbc_scores, bc_scores);
}

TEST_F(DynBetweennessGTest, runApproxBetweenness) {
    DorogovtsevMendesGenerator generator(100);
    Graph G1 = generator.generate();
    Graph G(G1, true, false);
    ApproxBetweenness bc(G, 0.1, 0.1);
    bc.run();
    DEBUG("Number of samples: ", bc.numberOfSamples());
    ApproxBetweenness bc1(G1, 0.1, 0.1);
    bc1.run();
    DEBUG("Number of samples: ", bc1.numberOfSamples());
}

TEST_F(DynBetweennessGTest, runDynVsStaticCaseInsertDirected){
    Aux::Random::setSeed(0, false);

    for(count n = 2; n <= 25; n++)
        for(count t = 0; t < 100; t++){
            auto g = ErdosRenyiGenerator(n, 0.3, true).generate();

            if(g.numberOfEdges() == g.numberOfNodes() * (g.numberOfNodes() - 1))continue;

            auto [u, v] = getNonAdjacentNodes(g);
            auto ibet = DynBetweenness(g);

            ibet.run();

            auto ne = GraphEvent(GraphEvent::EDGE_ADDITION, u, v);
            g.addEdge(u, v);
            ibet.update(ne);

            EXPECT_TRUE(g.hasEdge(u, v));
            auto brandes = Betweenness(g);
            brandes.run();
            g.forNodes([&](node w){
                EXPECT_NEAR(brandes.score(w), ibet.score(w), 1e-8);
            });
        }
}

TEST_F(DynBetweennessGTest, runDynVsStaticCaseInsertUndirected){
    Aux::Random::setSeed(0, false);

    for(count n = 2; n <= 25; n++)
        for(count t = 0; t < 100; t++){
            auto g = ErdosRenyiGenerator(n, 0.3, false).generate();
            if(g.numberOfEdges() == g.numberOfNodes() * (g.numberOfNodes() - 1) / 2)continue;

            auto [u, v] = getNonAdjacentNodes(g);
            auto ibet = DynBetweenness(g);

            ibet.run();

            auto ne = GraphEvent(GraphEvent::EDGE_ADDITION, u, v);
            g.addEdge(u, v);
            ibet.update(ne);
            EXPECT_TRUE(g.hasEdge(u, v));

            auto brandes = Betweenness(g);
            brandes.run();
            g.forNodes([&](node w){
                EXPECT_NEAR(brandes.score(w), ibet.score(w), 1e-8);
            });
        }
}

} /* namespace NetworKit */
