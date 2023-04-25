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
#include <networkit/generators/HyperbolicGenerator.hpp>

namespace NetworKit {

class ApproxElectricalClosenessGTest : public testing::Test {};

TEST_F(ApproxElectricalClosenessGTest, testApproxElectricalCloseness) {
    const double eps = 0.1;
    const count n_gen = 75;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);
        auto G = HyperbolicGenerator(n_gen, 10, 3).generate();
        G = ConnectedComponents::extractLargestConnectedComponent(G, true);
        count n = G.numberOfNodes();

        // Create a biconnected component with size 2.
        G.addNodes(2);
        G.addEdge(n - 1, n);
        G.addEdge(n, n + 1);

        ApproxElectricalCloseness apx(G);
        apx.run();
        const auto diag = apx.getDiagonal();
        const auto gt = apx.computeExactDiagonal(1e-12);
        G.forNodes([&](node u) { EXPECT_NEAR(diag[u], gt[u], eps); });
        EXPECT_EQ(apx.scores().size(), G.numberOfNodes());
    }
}

TEST_F(ApproxElectricalClosenessGTest, testDynApproxElectricalCloseness_run) {
    const double eps = 0.1;
    const count n_gen = 75;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);
        auto G = HyperbolicGenerator(n_gen, 10, 3).generate();
        G = ConnectedComponents::extractLargestConnectedComponent(G, true);
        count n = G.numberOfNodes();

        // Create a biconnected component with size 2.
        G.addNodes(2);
        G.addEdge(n - 1, n);
        G.addEdge(n, n + 1);

        ApproxElectricalCloseness apx(G);
        DynApproxElectricalCloseness dapx(G);
        apx.run();
        dapx.run();
        const auto diag = apx.getDiagonal();
        const auto ddiag = dapx.getDiagonal();
        G.forNodes([&](node u) { EXPECT_NEAR(diag[u], ddiag[u], eps); });
        EXPECT_EQ(dapx.scores().size(), G.numberOfNodes());
    }
}

TEST_F(ApproxElectricalClosenessGTest, testDynApproxElectricalCloseness_batchEdgeAddition) {
    const double eps = 0.1;
    count n_gen;
    for (int seed : {1, 2, 3}) {
        n_gen = 75;
        Aux::Random::setSeed(seed, true);
        auto G = HyperbolicGenerator(n_gen, 6, 3).generate();
        G = ConnectedComponents::extractLargestConnectedComponent(G, true);
        count n = G.numberOfNodes();

        // Create a biconnected component with size 2.
        G.addNodes(2);
        G.addEdge(n - 1, n);
        G.addEdge(n, n + 1);

        DynApproxElectricalCloseness dapx(G);
        dapx.run();

        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> distN(0, n + 1);

        node a, b;

        // add 10 random edges
        std::vector<GraphEvent> batch(10);
        for (count i = 0; i < 10; i++) {
            do {
                a = distN(rng);
                b = distN(rng);
            } while (G.hasEdge(a, b) || a == b);

            G.addEdge(a, b);

            batch[i].type = GraphEvent::EDGE_ADDITION;
            batch[i].u = a;
            batch[i].v = b;
        }

        dapx.updateBatch(batch);
        ApproxElectricalCloseness apx(G);
        apx.run();

        const auto diag = apx.getDiagonal();
        const auto ddiag = dapx.getDiagonal();
        G.forNodes([&](node u) { EXPECT_NEAR(diag[u], ddiag[u], eps); });
        EXPECT_EQ(dapx.scores().size(), G.numberOfNodes());
    }
}

TEST_F(ApproxElectricalClosenessGTest, testDynApproxElectricalCloseness_copy) {
    const double eps = 0.1;
    const count n_gen = 75;
    auto G = HyperbolicGenerator(n_gen, 10, 3).generate();
    G = ConnectedComponents::extractLargestConnectedComponent(G, true);
    count n = G.numberOfNodes();

    // Create a biconnected component with size 2.
    G.addNodes(2);
    G.addEdge(n - 1, n);
    G.addEdge(n, n + 1);

    DynApproxElectricalCloseness dapx(G);
    dapx.run();
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> distN(0, n + 1);

    node a, b;

    // add 10 random edges
    std::vector<GraphEvent> batch(10);
    for (count i = 0; i < 10; i++) {
        do {
            a = distN(rng);
            b = distN(rng);
        } while (G.hasEdge(a, b) || a == b);

        G.addEdge(a, b);

        batch[i].type = GraphEvent::EDGE_ADDITION;
        batch[i].u = a;
        batch[i].v = b;
    }

    DynApproxElectricalCloseness dapx2 = dapx;

    dapx2.updateBatch(batch);
    ApproxElectricalCloseness apx(G);
    apx.run();

    const auto diag = apx.getDiagonal();
    const auto ddiag2 = dapx2.getDiagonal();
    G.forNodes([&](node u) { EXPECT_NEAR(diag[u], ddiag2[u], eps); });
    EXPECT_EQ(dapx2.scores().size(), G.numberOfNodes());
}

} /* namespace NetworKit */