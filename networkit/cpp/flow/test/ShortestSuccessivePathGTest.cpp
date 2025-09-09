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

    template <typename CapacityAttr>
    static void addCostAndCapacity(Graph &G, CapacityAttr &capacities, node u, node v, double cost,
                                   double capacity) {
        G.addEdge(u, v, cost);
        capacities.set(G.edgeId(u, v), capacity);
    }
};

TEST_F(ShortestSuccessivePathGTest, testConstructorThrowsForUndirected) {
    Graph G(2, /*weighted*/ true, /*directed*/ false);
    try {
        MinFlowShortestSuccessivePath test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for undirected graph";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "MinFlowShortestSuccessivePath: Graph must be directed.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(ShortestSuccessivePathGTest, testConstructorThrowsForUnweighted) {
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

TEST_F(ShortestSuccessivePathGTest, testConstructorThrowsForUnindexedEdges) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    try {
        MinFlowShortestSuccessivePath test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for unindexed edges.";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "MinFlowShortestSuccessivePath: Graph edges must be indexed.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(ShortestSuccessivePathGTest, testConstructorThrowsForMissingEdgeAttribute) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    G.attachNodeDoubleAttribute(supplyName);

    try {
        MinFlowShortestSuccessivePath test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for missing edge attribute.";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), ("MinFlowShortestSuccessivePath: Provided edge attribute '"
                                + capacityName + "' not found.")
                                   .c_str());
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(ShortestSuccessivePathGTest, testConstructorThrowsForMissingNodeAttribute) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    G.attachEdgeDoubleAttribute(capacityName);

    try {
        MinFlowShortestSuccessivePath test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for missing node attribute.";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), ("MinFlowShortestSuccessivePath: Provided node attribute '"
                                + supplyName + "' not found.")
                                   .c_str());
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(ShortestSuccessivePathGTest, testConstructorThrowsForNegativeCapacity) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    G.addEdge(0, 1);
    auto eid = G.edgeId(0, 1);
    capacities.set(eid, -5.0); // invalid negative capacity
    supply.set(0, 0.0);
    supply.set(1, 0.0);

    try {
        MinFlowShortestSuccessivePath alg(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for negative capacity..";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "MinFlowShortestSuccessivePath: Capacities must be non-negative.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(ShortestSuccessivePathGTest, testConstructorThrowsForUnbalancedSupplies) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    G.addEdge(0, 1);
    const auto eid = G.edgeId(0, 1);
    capacities.set(eid, 5.0);
    supply.set(0, -5.0);
    supply.set(1, 4.0);

    try {
        MinFlowShortestSuccessivePath alg(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for unbalanced supplies/demands in nodes.";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "MinFlowShortestSuccessivePath: Sum of node supplies and demands "
                               "does not add up to zero.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(ShortestSuccessivePathGTest, testConstructorSucceedsWhenInputValid) {
    Graph G(3, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);
    for (node u = 0; u < 3; ++u) {
        supply.set(u, 0.0);
    }

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

TEST_F(ShortestSuccessivePathGTest, testRunThrowsOnNegativeCostCycle) {
    // 3-node directed graph with a negative cycle: 0->1->2->0, total cost = 1 + 1 - 5 = -3
    Graph G(3, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    // Add the negative cycle edges with positive capacities
    G.addEdge(0, 1, /*cost*/ +1.0);
    G.addEdge(1, 2, /*cost*/ +1.0);
    G.addEdge(2, 0, /*cost*/ -5.0); // completes a negative cycle

    for (node u = 0; u < 3; ++u) {
        supply.set(u, 0.0);
    }
    capacities.set(G.edgeId(0, 1), 3.0);
    capacities.set(G.edgeId(1, 2), 3.0);
    capacities.set(G.edgeId(2, 0), 3.0);

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

TEST_F(ShortestSuccessivePathGTest, testSimple5NodeGraph) {
    Graph G(5, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);

    auto supply = G.attachNodeDoubleAttribute(supplyName);
    for (node u = 0; u < 5; ++u) {
        supply.set(u, 0.0);
    }
    supply.set(0, +5.0);
    supply.set(4, -5.0);

    G.addEdge(0, 1, /*cost*/ 1.0);
    G.addEdge(1, 4, /*cost*/ 1.0);
    G.addEdge(0, 2, /*cost*/ 2.0);
    G.addEdge(2, 4, /*cost*/ 1.0);

    capacities.set(G.edgeId(0, 1), 3.0);
    capacities.set(G.edgeId(1, 4), 3.0);
    capacities.set(G.edgeId(0, 2), 3.0);
    capacities.set(G.edgeId(2, 4), 3.0);

    MinFlowShortestSuccessivePath solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 12.0);
}

TEST_F(ShortestSuccessivePathGTest, testComplexMultiSourceNegativeCost) {
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

    addCostAndCapacity(G, capacities, 0, 2, /*cost*/ 2.0, /*capacity*/ 4.0);
    addCostAndCapacity(G, capacities, 2, 5, /*cost*/ 1.0, /*capacity*/ 4.0);
    addCostAndCapacity(G, capacities, 0, 3, /*cost*/ 1.0, /*capacity*/ 2.0);
    addCostAndCapacity(G, capacities, 3, 5, /*cost*/ 3.0, /*capacity*/ 2.0);
    addCostAndCapacity(G, capacities, 1, 3, /*cost*/ -1.0, /*capacity*/ 3.0);
    addCostAndCapacity(G, capacities, 3, 4, /*cost*/ 1.0, /*capacity*/ 3.0);
    addCostAndCapacity(G, capacities, 4, 5, /*cost*/ 2.0, /*capacity*/ 3.0);
    addCostAndCapacity(G, capacities, 1, 2, /*cost*/ 3.0, /*capacity*/ 3.0);

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

    addCostAndCapacity(G, capacities, 0, 1, /*cost*/ 0.0, /*capacity*/ 4.0);
    addCostAndCapacity(G, capacities, 0, 2, /*cost*/ 0.0, /*capacity*/ 4.0);
    addCostAndCapacity(G, capacities, 4, 3, /*cost*/ 0.0, /*capacity*/ 4.0);
    addCostAndCapacity(G, capacities, 5, 3, /*cost*/ 0.0, /*capacity*/ 4.0);
    addCostAndCapacity(G, capacities, 0, 3, /*cost*/ 0.0, /*capacity*/ 0.0);

    MinFlowShortestSuccessivePath solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 0.0);
}

// https://lpsolve.sourceforge.net/5.5/DIMACS_mcf.htm
TEST_F(ShortestSuccessivePathGTest, testSimpleDimacsProblem) {
    Graph G(4, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    supply.set(0, 4.0);
    supply.set(1, 0.0);
    supply.set(2, 0.0);
    supply.set(3, -4.0);

    addCostAndCapacity(G, capacities, 0, 2, /*cost*/ 2.0, /*capacity*/ 2.0);
    addCostAndCapacity(G, capacities, 0, 1, /*cost*/ 2.0, /*capacity*/ 4.0);
    addCostAndCapacity(G, capacities, 1, 2, /*cost*/ 1.0, /*capacity*/ 2.0);
    addCostAndCapacity(G, capacities, 1, 3, /*cost*/ 3.0, /*capacity*/ 3.0);
    addCostAndCapacity(G, capacities, 2, 3, /*cost*/ 1.0, /*capacity*/ 5.0);

    MinFlowShortestSuccessivePath solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 14.0);
}

TEST_F(ShortestSuccessivePathGTest, testSwapNeededByBackwardResiduals) {
    Graph G(6, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    // Supply / demand: send 2 units from 0 to 5
    supply.set(0, 2.0);
    supply.set(1, 0.0);
    supply.set(2, 0.0);
    supply.set(3, 0.0);
    supply.set(4, 0.0);
    supply.set(5, -2.0);

    addCostAndCapacity(G, capacities, 0, 1, /*cost*/ 0.0, /*capacity*/ 1.0);
    addCostAndCapacity(G, capacities, 0, 2, /*cost*/ 0.0, /*capacity*/ 1.0);

    addCostAndCapacity(G, capacities, 1, 3, /*cost*/ 1.0, /*capacity*/ 1.0);
    addCostAndCapacity(G, capacities, 1, 4, /*cost*/ 100.0, /*capacity*/ 1.0);
    addCostAndCapacity(G, capacities, 2, 3, /*cost*/ 0.0, /*capacity*/ 1.0);
    addCostAndCapacity(G, capacities, 2, 4, /*cost*/ 0.0, /*capacity*/ 1.0);

    addCostAndCapacity(G, capacities, 3, 5, /*cost*/ 0.0, /*capacity*/ 1.0);
    addCostAndCapacity(G, capacities, 4, 5, /*cost*/ 0.0, /*capacity*/ 1.0);

    MinFlowShortestSuccessivePath solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 1.0);
}

TEST_F(ShortestSuccessivePathGTest, testMultipleSuppliesAndDemands) {
    Graph G(6, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    supply.set(0, 2.0);
    supply.set(1, -1.0);
    supply.set(2, -1.0);
    supply.set(3, 4.0);
    supply.set(4, -3.0);
    supply.set(5, -1.0);

    addCostAndCapacity(G, capacities, 0, 1, /*cost*/ 2.0, /*capacity*/ 1.0);
    addCostAndCapacity(G, capacities, 0, 2, /*cost*/ 3.0, /*capacity*/ 2.0);

    addCostAndCapacity(G, capacities, 0, 3, /*cost*/ 0.0, /*capacity*/ 1.0);
    addCostAndCapacity(G, capacities, 1, 3, /*cost*/ 0.0, /*capacity*/ 1.0);
    addCostAndCapacity(G, capacities, 2, 5, /*cost*/ 0.0, /*capacity*/ 1.0);

    addCostAndCapacity(G, capacities, 3, 4, /*cost*/ 3.0, /*capacity*/ 10.0);
    addCostAndCapacity(G, capacities, 3, 5, /*cost*/ 1.0, /*capacity*/ 1.0);

    MinFlowShortestSuccessivePath solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 15.0);
}

TEST_F(ShortestSuccessivePathGTest, testMultipleSuppliesAndDemandsWithDetour) {
    Graph G(7, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    supply.set(0, 21.0);
    supply.set(1, -1.0);
    supply.set(2, -1.0);
    supply.set(3, 21.0);
    supply.set(4, -3.0);
    supply.set(5, 0.0);
    supply.set(6, -37.0);

    addCostAndCapacity(G, capacities, 0, 3, /*cost*/ 0.0, /*capacity*/ 1.0);
    addCostAndCapacity(G, capacities, 0, 5, /*cost*/ 2.0, /*capacity*/ 10.0);
    addCostAndCapacity(G, capacities, 0, 1, /*cost*/ 2.0, /*capacity*/ 1.0);
    addCostAndCapacity(G, capacities, 0, 2, /*cost*/ 3.0, /*capacity*/ 10.0);

    addCostAndCapacity(G, capacities, 2, 5, /*cost*/ 2.0, /*capacity*/ 12.0);
    addCostAndCapacity(G, capacities, 1, 3, /*cost*/ 0.0, /*capacity*/ 1.0);

    addCostAndCapacity(G, capacities, 3, 4, /*cost*/ 3.0, /*capacity*/ 50.0);
    addCostAndCapacity(G, capacities, 3, 5, /*cost*/ 1.0, /*capacity*/ 1.0);

    addCostAndCapacity(G, capacities, 5, 6, /*cost*/ 4.0, /*capacity*/ 45.0);
    addCostAndCapacity(G, capacities, 4, 6, /*cost*/ 4.0, /*capacity*/ 20.0);

    MinFlowShortestSuccessivePath solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 277.0);
}

} // namespace NetworKit
