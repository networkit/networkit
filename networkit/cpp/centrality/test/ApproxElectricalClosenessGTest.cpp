/*
 * ApproxElectricalClosenessGTest.cpp
 *
 *  Created on: 18.04.2023
 *      Author: Lukas Berner <Lukas.Berner@hu-berlin.de>
 */

#include <gtest/gtest.h>

#include <networkit/centrality/ApproxElectricalCloseness.hpp>
#include <networkit/centrality/DynApproxElectricalCloseness.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/generators/HyperbolicGenerator.hpp>
#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

class ApproxElectricalClosenessGTest
    : public testing::TestWithParam<std::tuple<int, std::string, double>> {
protected:
    const count n_gen = 75;
    const double er_prob = 0.15;

    Graph generate(std::string generator) {
        Graph G;
        if (generator == "ER")
            G = ErdosRenyiGenerator(n_gen, er_prob).generate();
        else if (generator == "Hyperbolic")
            G = HyperbolicGenerator(n_gen).generate();
        else
            throw std::logic_error("unkown generator");
        G = ConnectedComponents::extractLargestConnectedComponent(G, true);
        const count n = G.numberOfNodes();

        // Create a biconnected component with size 2.
        G.addNodes(2);
        G.addEdge(n - 1, n);
        G.addEdge(n, n + 1);

        return G;
    };

    GraphEvent add_random_edge(Graph &G) {
        node a, b;

        do {
            a = Aux::Random::integer(0, G.numberOfNodes() - 1);
            b = Aux::Random::integer(0, G.numberOfNodes() - 1);
        } while (G.hasEdge(a, b) || a == b);

        G.addEdge(a, b);
        return GraphEvent(GraphEvent::EDGE_ADDITION, a, b);
    };

    GraphEvent remove_random_edge(Graph &G, const node pivot = none) {
        node a, b;

        do {
            a = Aux::Random::integer(0, G.numberOfNodes() - 1);
            b = Aux::Random::integer(0, G.numberOfNodes() - 1);
            if (a > b)
                std::swap(a, b);
        } while (!G.hasEdge(a, b) || a == pivot || b == pivot);

        G.removeEdge(a, b);
        return GraphEvent(GraphEvent::EDGE_REMOVAL, a, b);
    };
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, ApproxElectricalClosenessGTest,
                         testing::Combine(testing::Values(1, 2, 3),            // seed
                                          testing::Values("ER", "Hyperbolic"), // generator
                                          testing::Values(0.1, 0.3, 0.5)       // eps
                                          ));

TEST_P(ApproxElectricalClosenessGTest, testApproxElectricalCloseness) {
    auto [seed, generator, eps] = GetParam();
    Aux::Random::setSeed(seed, true);
    auto G = generate(generator);

    ApproxElectricalCloseness apx(G, eps);
    apx.run();
    const auto diag = apx.getDiagonal();
    const auto gt = apx.computeExactDiagonal(1e-12);
    G.forNodes([&eps = eps, &gt, &diag](node u) { EXPECT_NEAR(diag[u], gt[u], eps); });
    EXPECT_EQ(apx.scores().size(), G.numberOfNodes());
}

TEST_P(ApproxElectricalClosenessGTest, testDynApproxElectricalCloseness_run) {
    auto [seed, generator, eps] = GetParam();
    Aux::Random::setSeed(seed, true);
    auto G = generate(generator);

    ApproxElectricalCloseness apx(G, eps);
    DynApproxElectricalCloseness dapx(G, eps);
    apx.run();
    dapx.run();
    const auto diag = apx.getDiagonal();
    const auto ddiag = dapx.getDiagonal();
    G.forNodes([&eps = eps, &ddiag, &diag](node u) {
        EXPECT_NEAR(diag[u], ddiag[u], 2 * eps);
    }); // 2 * eps because both have max abs error of eps
    EXPECT_EQ(dapx.scores().size(), G.numberOfNodes());
}

TEST_P(ApproxElectricalClosenessGTest, testDynApproxElectricalCloseness_EdgeAddition) {
    auto [seed, generator, eps] = GetParam();
    Aux::Random::setSeed(seed, true);
    auto G = generate(generator);

    DynApproxElectricalCloseness dapx(G, eps);
    dapx.run();

    for (int i = 0; i < 10; i++) {
        auto event = add_random_edge(G);
        dapx.update(event);
    }

    ApproxElectricalCloseness apx(G, eps);
    apx.run();

    const auto diag = apx.getDiagonal();
    const auto ddiag = dapx.getDiagonal();
    G.forNodes([&eps = eps, &ddiag, &diag](node u) { EXPECT_NEAR(diag[u], ddiag[u], 2 * eps); });
    EXPECT_EQ(dapx.scores().size(), G.numberOfNodes());
}

TEST_P(ApproxElectricalClosenessGTest, testDynApproxElectricalCloseness_EdgeDeletion) {
    omp_set_num_threads(1);
    auto [seed, generator, eps] = GetParam();
    Aux::Random::setSeed(seed, true);
    auto G = generate(generator);

    // make sure that removing an edge does not disconnect the graph
    auto star_node = G.addNode();
    G.forNodes([&](node u) {
        if (u != star_node)
            G.addEdge(star_node, u);
    });

    DynApproxElectricalCloseness dapx(G, eps, 0.3, star_node);
    dapx.run();

    auto event = remove_random_edge(G, star_node);
    dapx.update(event);

    ApproxElectricalCloseness apx(G, eps);
    apx.run();

    const auto diag = apx.getDiagonal();
    const auto ddiag = dapx.getDiagonal();
    G.forNodes([&eps = eps, &ddiag, &diag](node u) { EXPECT_NEAR(diag[u], ddiag[u], 2 * eps); });
    EXPECT_EQ(dapx.scores().size(), G.numberOfNodes());
}

TEST_P(ApproxElectricalClosenessGTest, testDynApproxElectricalCloseness_copy) {
    auto [seed, generator, eps] = GetParam();
    Aux::Random::setSeed(seed, true);
    auto G = generate(generator);

    DynApproxElectricalCloseness dapx(G, eps);
    dapx.run();

    auto event = add_random_edge(G);

    // make copy
    DynApproxElectricalCloseness dapx2 = dapx;

    dapx2.update(event);
    ApproxElectricalCloseness apx(G, eps);
    apx.run();

    const auto diag = apx.getDiagonal();
    const auto ddiag2 = dapx2.getDiagonal();
    G.forNodes([&eps = eps, &ddiag2, &diag](node u) { EXPECT_NEAR(diag[u], ddiag2[u], 2 * eps); });
    EXPECT_EQ(dapx2.scores().size(), G.numberOfNodes());
}

} /* namespace NetworKit */