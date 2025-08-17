/*  ShortestSuccessivePath.hpp
 *
 *	Created on: 05.08.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <gtest/gtest.h>
#include <networkit/flow/ShortestSuccessivePath.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class ShortestSuccessivePathGTest : public ::testing::Test {
protected:
    const std::string capacityName = "capacity";
    const std::string supplyName = "supply";
};

TEST_F(ShortestSuccessivePathGTest, constructorThrowsForUndirected) {
    Graph G(2, /*weighted*/ true, /*directed*/ false);
    try {
        MinFlowShortestSuccessivePath test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for undirected graph";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "MinFlowShortestSuccessivePath: Graph must be directed");
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(ShortestSuccessivePathGTest, constructorThrowsForUnweighted) {
    Graph G(2, /*weighted*/ false, /*directed*/ true);
    try {
        MinFlowShortestSuccessivePath test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for unweighted graph";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "MinFlowShortestSuccessivePath: Graph must be weighted.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(ShortestSuccessivePathGTest, constructorThrowsForUnindexedEdges) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    try {
        MinFlowShortestSuccessivePath test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for unindexed edges";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "MinFlowShortestSuccessivePath: Graph edges must be indexed");
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(ShortestSuccessivePathGTest, constructorThrowsForMissingEdgeAttribute) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    G.attachNodeDoubleAttribute(supplyName);

    try {
        MinFlowShortestSuccessivePath test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for missing edge attribute";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), ("MinFlowShortestSuccessivePath: Provided edge attribute '"
                                + capacityName + "' not found")
                                   .c_str());
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(ShortestSuccessivePathGTest, constructorThrowsForMissingNodeAttribute) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    G.attachEdgeDoubleAttribute(capacityName);

    try {
        MinFlowShortestSuccessivePath test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for missing node attribute";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), ("MinFlowShortestSuccessivePath: Provided node attribute '"
                                + supplyName + "' not found")
                                   .c_str());
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(ShortestSuccessivePathGTest, constructorThrowsForNegativeCapacity) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();

    auto capacities   = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    G.addEdge(0, 1);
    auto eid = G.edgeId(0, 1);
    capacities.set(eid, -5.0);  // invalid negative capacity

    try {
        MinFlowShortestSuccessivePath alg(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for negative capacity";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(),
                     "MinFlowShortestSuccessivePath: Capacities must be non-negative");
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(ShortestSuccessivePathGTest, constructorSucceedsWhenValid) {
    Graph G(3, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    EXPECT_NO_THROW({ MinFlowShortestSuccessivePath test(G, capacityName, supplyName); });
}

TEST_F(ShortestSuccessivePathGTest, testZeroNodesGraph) {
    Graph G(0, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    MinFlowShortestSuccessivePath test(G, capacityName, supplyName);
    test.run();
    EXPECT_DOUBLE_EQ(test.getTotalCost(), 0.0);
}

TEST_F(ShortestSuccessivePathGTest, runThrowsOnNegativeCostCycle) {
    // 3-node directed graph with a negative cycle: 0->1->2->0, total cost = 1 + 1 - 5 = -3
    Graph G(3, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply     = G.attachNodeDoubleAttribute(supplyName);

    // Add the negative cycle edges with positive capacities
    G.addEdge(0, 1, /*cost=*/+1.0);
    G.addEdge(1, 2, /*cost=*/+1.0);
    G.addEdge(2, 0, /*cost=*/-5.0); // completes a negative cycle

    for (node u = 0; u < 3; ++u) {
        supply.set(u, 0.0);
    }
    capacities.set(G.edgeId(0,1), 3.0);
    capacities.set(G.edgeId(1,2), 3.0);
    capacities.set(G.edgeId(2,0), 3.0);

    MinFlowShortestSuccessivePath solver(G, capacityName, supplyName);

    try {
        solver.run();
        FAIL() << "Expected std::runtime_error due to negative-cost cycle";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(),
            "MinFlowShortestSuccessivePath: negative-cost cycle in residual graph");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(ShortestSuccessivePathGTest, Simple5NodeGraph) {
    Graph G(5, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);

    auto supply = G.attachNodeDoubleAttribute(supplyName);
    for (node u = 0; u < 5; ++u) {
        supply.set(u, 0.0);
    }
    supply.set(0, +5.0);
    supply.set(4, -5.0);

    G.addEdge(0, 1, /*cost=*/1.0);
    G.addEdge(1, 4, /*cost=*/1.0);
    G.addEdge(0, 2, /*cost=*/2.0);
    G.addEdge(2, 4, /*cost=*/1.0);

    capacities.set(G.edgeId(0,1), 3.0);
    capacities.set(G.edgeId(1,4), 3.0);
    capacities.set(G.edgeId(0,2), 3.0);
    capacities.set(G.edgeId(2,4), 3.0);

    MinFlowShortestSuccessivePath solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 12.0);
}

TEST_F(ShortestSuccessivePathGTest, ComplexMultiSourceNegativeCost) {
    Graph G(6, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    for (node u = 0; u < 6; ++u) {
        supply.set(u, 0.0);
    }
    supply.set(0, 4.0);
    supply.set(1, 3.0);
    supply.set(5, -7.0);

    auto addCostAndCapacity = [&](node u, node v, double cost, double capacity) {
        G.addEdge(u, v, cost);
        capacities.set(G.edgeId(u, v), capacity);
    };

    addCostAndCapacity(0, 2, /*cost=*/2.0, /*capacity=*/4.0);
    addCostAndCapacity(2, 5, /*cost=*/1.0, /*capacity=*/4.0);
    addCostAndCapacity(0, 3, /*cost=*/1.0, /*capacity=*/2.0);
    addCostAndCapacity(3, 5, /*cost=*/3.0, /*capacity=*/2.0);
    addCostAndCapacity(1, 3, /*cost=*/-1.0,/*capacity=*/3.0);
    addCostAndCapacity(3, 4, /*cost=*/1.0, /*capacity=*/3.0);
    addCostAndCapacity(4, 5, /*cost=*/2.0, /*capacity=*/3.0);
    addCostAndCapacity(1, 2, /*cost=*/3.0, /*capacity=*/3.0);

    MinFlowShortestSuccessivePath solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 18.0);
}

TEST_F(ShortestSuccessivePathGTest, testZeroWeights) {
    Graph G(6, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    supply.set(0, 4.0);
    supply.set(1, -2.0);
    supply.set(2, -2.0);
    supply.set(3, -4.0);
    supply.set(4, +2.0);
    supply.set(5, +2.0);

    auto addZeroWeightAndCapacity = [&](node u, node v, double capacity) {
        G.addEdge(u, v, 0.0);
        capacities.set(G.edgeId(u, v), capacity);
    };

    addZeroWeightAndCapacity(0, 1, /*capacity=*/4.0);
    addZeroWeightAndCapacity(0, 2,/*capacity=*/4.0);
    addZeroWeightAndCapacity(4, 3, /*capacity=*/4.0);
    addZeroWeightAndCapacity(5, 3,/*capacity=*/4.0);
    addZeroWeightAndCapacity(0, 3, /*capacity=*/0.0);

    MinFlowShortestSuccessivePath solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 0.0);
}


} // namespace NetworKit
