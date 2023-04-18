/*
 * ApproxElectricalClosenessGTest.cpp
 *
 *  Created on: 18.04.2023
 *      Author: Lukas Berner <Lukas.Berner@hu-berlin.de>
 */

#include <gtest/gtest.h>

#include <networkit/centrality/ApproxElectricalCloseness.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/generators/HyperbolicGenerator.hpp>

namespace NetworKit {

class ApproxElectricalClosenessGTest : public testing::Test {};

TEST_F(ApproxElectricalClosenessGTest, testApproxElectricalCloseness) {
    const double eps = 0.1;
    const count n = 75;
    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, true);
        auto G = HyperbolicGenerator(n, 10, 3).generate();
        G = ConnectedComponents::extractLargestConnectedComponent(G, true);

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

} /* namespace NetworKit */