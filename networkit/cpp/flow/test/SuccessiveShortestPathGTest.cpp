/*  ShortestSuccessivePath.hpp
 *
 *	Created on: 05.08.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <networkit/flow/SuccessiveShortestPath.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
class SuccessiveShortestPathGTest : public ::testing::Test {
protected:
    const std::string capacityName = "capacity";
    const std::string supplyName = "supply";

    template <typename CapacityAttr>
    static void addCostAndCapacity(Graph &G, CapacityAttr &capacities, node u, node v, double cost,
                                   double capacity) {
        G.addEdge(u, v, cost);
        capacities.set(G.edgeId(u, v), capacity);
    }

    // For every node: (outgoing flow − incoming flow) must equal its supply/demand.
    void checkFlowConservation(Graph &G, auto flow) {
        std::vector<double> imbalance(G.numberOfNodes(), 0.0);
        const auto supply = G.getNodeDoubleAttribute(supplyName);
        G.forEdges([&](node u, node v, edgeid eid) {
            const double f = flow.get(eid);
            imbalance[u] += f; // outgoing flow
            imbalance[v] -= f; // incoming flow
        });
        G.forNodes([&](node u) { EXPECT_DOUBLE_EQ(imbalance[u], supply.get(u)); });
    }
};

TEST_F(SuccessiveShortestPathGTest, testConstructorThrowsForUndirected) {
    Graph G(2, /*weighted*/ true, /*directed*/ false);
    try {
        SuccessiveShortestPathMinCostFlow test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for undirected graph";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "SuccessiveShortestPathMinCostFlow: Graph must be directed.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(SuccessiveShortestPathGTest, testConstructorThrowsForUnweighted) {
    Graph G(2, /*weighted*/ false, /*directed*/ true);
    try {
        SuccessiveShortestPathMinCostFlow test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for unweighted graph";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "SuccessiveShortestPathMinCostFlow: Graph must be weighted.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(SuccessiveShortestPathGTest, testConstructorThrowsForUnindexedEdges) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    try {
        SuccessiveShortestPathMinCostFlow test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for unindexed edges.";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "SuccessiveShortestPathMinCostFlow: Graph edges must be indexed.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(SuccessiveShortestPathGTest, testConstructorThrowsForMissingEdgeAttribute) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    G.attachNodeDoubleAttribute(supplyName);

    try {
        SuccessiveShortestPathMinCostFlow test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for missing edge attribute.";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), ("SuccessiveShortestPathMinCostFlow: Provided edge attribute '"
                                + capacityName + "' not found.")
                                   .c_str());
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(SuccessiveShortestPathGTest, testConstructorThrowsForMissingNodeAttribute) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    G.attachEdgeDoubleAttribute(capacityName);

    try {
        SuccessiveShortestPathMinCostFlow test(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for missing node attribute.";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), ("SuccessiveShortestPathMinCostFlow: Provided node attribute '"
                                + supplyName + "' not found.")
                                   .c_str());
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(SuccessiveShortestPathGTest, testConstructorThrowsForNegativeCapacity) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    G.addEdge(0, 1);
    auto eid = G.edgeId(0, 1);
    capacities.set(eid, -5.0);
    supply.set(0, 0.0);
    supply.set(1, 0.0);

    try {
        SuccessiveShortestPathMinCostFlow alg(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for negative capacity..";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(),
                     "SuccessiveShortestPathMinCostFlow: Capacities must be non-negative.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(SuccessiveShortestPathGTest, testConstructorThrowsForUnbalancedSupplies) {
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    G.addEdge(0, 1);
    const edgeid eid = G.edgeId(0, 1);
    capacities.set(eid, 5.0);
    supply.set(0, -5.0);
    supply.set(1, 4.0);

    try {
        SuccessiveShortestPathMinCostFlow alg(G, capacityName, supplyName);
        FAIL() << "Expected std::runtime_error for unbalanced supplies/demands in nodes.";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(),
                     "SuccessiveShortestPathMinCostFlow: Sum of node supplies and demands "
                     "does not add up to zero.");
    } catch (...) {
        FAIL() << "Expected std::runtime_error";
    }
}

TEST_F(SuccessiveShortestPathGTest, testConstructorSucceedsWhenInputValid) {
    Graph G(3, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);
    for (node u = 0; u < 3; ++u) {
        supply.set(u, 0.0);
    }

    EXPECT_NO_THROW({ SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName); });
}

TEST_F(SuccessiveShortestPathGTest, testRunNotCallesGetTotalCost) {
    Graph G(3, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);
    for (node u = 0; u < 3; ++u) {
        supply.set(u, 0.0);
    }
    SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName);
    try {
        solver.getTotalCost();
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Error, run must be called first");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(SuccessiveShortestPathGTest, testRunNotCallesGetFlow) {
    Graph G(3, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);
    for (node u = 0; u < 3; ++u) {
        supply.set(u, 0.0);
    }
    SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName);
    try {
        solver.getFlow();
        FAIL() << "Expected std::runtime_error";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(), "Error, run must be called first");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(SuccessiveShortestPathGTest, testZeroNodesGraph) {
    Graph G(0, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName);
    solver.run();
    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 0.0);
    const auto flow = solver.getFlow();
    EXPECT_THAT(flow, ::testing::SizeIs(0));
}

TEST_F(SuccessiveShortestPathGTest, testRunThrowsOnNegativeCostCycle) {
    Graph G(3, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    for (node u = 0; u < 3; ++u) {
        supply.set(u, 0.0);
    }

    addCostAndCapacity(G, capacities, 0, 1, /*cost*/ 1.0, /*capacity*/ 3.0);
    addCostAndCapacity(G, capacities, 1, 2, /*cost*/ 1.0, /*capacity*/ 3.0);
    addCostAndCapacity(G, capacities, 2, 0, /*cost*/ -5.0, /*capacity*/ 3.0);

    SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName);

    try {
        solver.run();
        FAIL() << "Expected std::runtime_error due to negative-cost cycle";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(),
                     "SuccessiveShortestPathMinCostFlow: negative-cost cycle in residual graph");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(SuccessiveShortestPathGTest, testRunThrowsWhenSuppliesDemandsCanNotBeSatisfied) {
    Graph G(3, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    supply.set(0, 10.0);
    supply.set(1, -5.0);
    supply.set(2, -5.0);

    addCostAndCapacity(G, capacities, 0, 1, /*cost*/ 1.0, /*capacity*/ 3.0);
    addCostAndCapacity(G, capacities, 1, 2, /*cost*/ 1.0, /*capacity*/ 3.0);
    addCostAndCapacity(G, capacities, 2, 0, /*cost*/ 2.0, /*capacity*/ 3.0);

    SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName);

    try {
        solver.run();
        FAIL() << "Expected std::runtime_error because supplies/demands can not be satisfied by "
                  "network.";
    } catch (const std::runtime_error &e) {
        EXPECT_STREQ(e.what(),
                     "SuccessiveShortestPathMinCostFlow: unable to satisfy all supplies/demands");
    } catch (...) {
        FAIL() << "Expected std::runtime_error but got a different exception.";
    }
}

TEST_F(SuccessiveShortestPathGTest, testSimple5NodeGraph) {
    Graph G(5, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);

    auto supply = G.attachNodeDoubleAttribute(supplyName);
    for (node u = 0; u < 5; ++u) {
        supply.set(u, 0.0);
    }
    supply.set(0, +5.0);
    supply.set(4, -5.0);

    addCostAndCapacity(G, capacities, 0, 1, /*cost*/ 1.0, /*capacity*/ 3.0);
    addCostAndCapacity(G, capacities, 1, 4, /*cost*/ 1.0, /*capacity*/ 3.0);
    addCostAndCapacity(G, capacities, 0, 2, /*cost*/ 2.0, /*capacity*/ 3.0);
    addCostAndCapacity(G, capacities, 2, 4, /*cost*/ 1.0, /*capacity*/ 3.0);

    SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 12.0);

    const auto flow = solver.getFlow();
    EXPECT_DOUBLE_EQ(flow.get(G.edgeId(0, 1)), 3.0);
    EXPECT_DOUBLE_EQ(flow.get(G.edgeId(1, 4)), 3.0);
    EXPECT_DOUBLE_EQ(flow.get(G.edgeId(0, 2)), 2.0);
    EXPECT_DOUBLE_EQ(flow.get(G.edgeId(2, 4)), 2.0);

    checkFlowConservation(G, flow);
}

TEST_F(SuccessiveShortestPathGTest, testComplexMultiSourceNegativeCost) {
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

    SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 18.0);

    const auto flow = solver.getFlow();
    std::vector<double> expected(G.upperEdgeIdBound(), 0.0);
    expected[G.edgeId(0, 2)] = 4.0;
    expected[G.edgeId(2, 5)] = 4.0;

    expected[G.edgeId(1, 3)] = 3.0;
    expected[G.edgeId(3, 5)] = 2.0;
    expected[G.edgeId(3, 4)] = 1.0;
    expected[G.edgeId(4, 5)] = 1.0;

    G.forEdges([&](node, node, edgeid eid) { EXPECT_DOUBLE_EQ(flow.get(eid), expected[eid]); });

    checkFlowConservation(G, flow);
}

TEST_F(SuccessiveShortestPathGTest, testZeroCosts) {
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

    SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 0.0);

    const auto flow = solver.getFlow();

    std::vector<double> expected(G.upperEdgeIdBound(), 0.0);
    expected[G.edgeId(0, 1)] = 2.0;
    expected[G.edgeId(0, 2)] = 2.0;
    expected[G.edgeId(4, 3)] = 2.0;
    expected[G.edgeId(5, 3)] = 2.0;

    G.forEdges([&](node, node, edgeid eid) { EXPECT_DOUBLE_EQ(flow.get(eid), expected[eid]); });

    checkFlowConservation(G, flow);
}

// https://lpsolve.sourceforge.net/5.5/DIMACS_mcf.htm
TEST_F(SuccessiveShortestPathGTest, testSimpleDimacsProblem) {
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

    SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 14.0);

    const auto flow = solver.getFlow();

    std::vector<double> expected(G.upperEdgeIdBound(), 0.0);
    expected[G.edgeId(0, 2)] = 2.0;
    expected[G.edgeId(2, 3)] = 4.0;
    expected[G.edgeId(0, 1)] = 2.0;
    expected[G.edgeId(1, 2)] = 2.0;

    G.forEdges([&](node, node, edgeid eid) { EXPECT_DOUBLE_EQ(flow.get(eid), expected[eid]); });

    checkFlowConservation(G, flow);
}

TEST_F(SuccessiveShortestPathGTest, testSwapNeededByBackwardResiduals) {
    Graph G(6, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

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

    SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 1.0);

    const auto flow = solver.getFlow();

    std::vector<double> expected(G.upperEdgeIdBound(), 0.0);
    expected[G.edgeId(0, 1)] = 1.0;
    expected[G.edgeId(0, 2)] = 1.0;

    expected[G.edgeId(1, 3)] = 1.0;
    expected[G.edgeId(1, 4)] = 0.0;
    expected[G.edgeId(2, 3)] = 0.0;
    expected[G.edgeId(2, 4)] = 1.0;

    expected[G.edgeId(3, 5)] = 1.0;
    expected[G.edgeId(4, 5)] = 1.0;

    G.forEdges([&](node, node, edgeid eid) { EXPECT_DOUBLE_EQ(flow.get(eid), expected[eid]); });

    checkFlowConservation(G, flow);
}

TEST_F(SuccessiveShortestPathGTest, testMultipleSuppliesAndDemands) {
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

    SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 15.0);
    const auto flow = solver.getFlow();

    std::vector<double> expected(G.upperEdgeIdBound(), 0.0);
    expected[G.edgeId(0, 1)] = 1.0;
    expected[G.edgeId(0, 2)] = 1.0;

    expected[G.edgeId(0, 3)] = 0.0;
    expected[G.edgeId(1, 3)] = 0.0;
    expected[G.edgeId(2, 5)] = 0.0;

    expected[G.edgeId(3, 4)] = 3.0;
    expected[G.edgeId(3, 5)] = 1.0;

    G.forEdges([&](node, node, edgeid eid) { EXPECT_DOUBLE_EQ(flow.get(eid), expected[eid]); });

    checkFlowConservation(G, flow);
}

TEST_F(SuccessiveShortestPathGTest, testMultipleSuppliesAndDemandsWithDetour) {
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

    SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 277.0);

    const auto flow = solver.getFlow();

    std::vector<double> expected(G.upperEdgeIdBound(), 0.0);
    expected[G.edgeId(0, 1)] = 1.0;
    expected[G.edgeId(0, 2)] = 9.0;
    expected[G.edgeId(0, 3)] = 1.0;
    expected[G.edgeId(0, 5)] = 10.0;

    expected[G.edgeId(1, 3)] = 0.0;
    expected[G.edgeId(2, 5)] = 8.0;

    expected[G.edgeId(3, 4)] = 21.0;
    expected[G.edgeId(3, 5)] = 1.0;

    expected[G.edgeId(4, 6)] = 18.0;
    expected[G.edgeId(5, 6)] = 19.0;

    G.forEdges([&](node, node, edgeid eid) { EXPECT_DOUBLE_EQ(flow.get(eid), expected[eid]); });
    checkFlowConservation(G, flow);
}

// generated with netgen http://archive.dimacs.rutgers.edu/pub/netflow/generators/network/netgen/
// using these inputs:
// ./netgen <<EOF > min10x20.min
// 12345
// 1
// 10 3 3 20 1 10 20 1 1 10 100 1 10
// EOF
TEST_F(SuccessiveShortestPathGTest, testDimacsGeneratedMin10x20) {
    Graph G(10, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    supply.set(0, 6.0);
    supply.set(1, 1.0);
    supply.set(2, 13.0);
    supply.set(3, 0.0);
    supply.set(4, 0.0);
    supply.set(5, 0.0);
    supply.set(6, 0.0);
    supply.set(7, -2.0);
    supply.set(8, -8.0);
    supply.set(9, -10.0);

    addCostAndCapacity(G, capacities, 0, 5, /*cost*/ 7.0, /*capacity*/ 6.0);
    addCostAndCapacity(G, capacities, 0, 8, /*cost*/ 8.0, /*capacity*/ 8.0);
    addCostAndCapacity(G, capacities, 0, 2, /*cost*/ 3.0, /*capacity*/ 2.0);

    addCostAndCapacity(G, capacities, 5, 6, /*cost*/ 6.0, /*capacity*/ 6.0);
    addCostAndCapacity(G, capacities, 5, 9, /*cost*/ 9.0, /*capacity*/ 6.0);
    addCostAndCapacity(G, capacities, 5, 8, /*cost*/ 7.0, /*capacity*/ 6.0);
    addCostAndCapacity(G, capacities, 5, 2, /*cost*/ 2.0, /*capacity*/ 7.0);
    addCostAndCapacity(G, capacities, 5, 4, /*cost*/ 6.0, /*capacity*/ 10.0);

    addCostAndCapacity(G, capacities, 6, 7, /*cost*/ 4.0, /*capacity*/ 6.0);

    addCostAndCapacity(G, capacities, 1, 3, /*cost*/ 6.0, /*capacity*/ 1.0);
    addCostAndCapacity(G, capacities, 1, 4, /*cost*/ 3.0, /*capacity*/ 9.0);
    addCostAndCapacity(G, capacities, 1, 8, /*cost*/ 6.0, /*capacity*/ 5.0);

    addCostAndCapacity(G, capacities, 3, 7, /*cost*/ 2.0, /*capacity*/ 1.0);
    addCostAndCapacity(G, capacities, 3, 9, /*cost*/ 10.0, /*capacity*/ 1.0);

    addCostAndCapacity(G, capacities, 2, 4, /*cost*/ 10.0, /*capacity*/ 13.0);

    addCostAndCapacity(G, capacities, 4, 9, /*cost*/ 5.0, /*capacity*/ 13.0);
    addCostAndCapacity(G, capacities, 4, 8, /*cost*/ 1.0, /*capacity*/ 13.0);

    addCostAndCapacity(G, capacities, 7, 9, /*cost*/ 2.0, /*capacity*/ 10.0);
    addCostAndCapacity(G, capacities, 7, 5, /*cost*/ 4.0, /*capacity*/ 7.0);
    addCostAndCapacity(G, capacities, 7, 6, /*cost*/ 7.0, /*capacity*/ 3.0);

    SuccessiveShortestPathMinCostFlow solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 248.0);

    const auto flow = solver.getFlow();
    std::vector<double> expected(G.upperEdgeIdBound(), 0.0);

    expected[G.edgeId(0, 5)] = 1.0;
    expected[G.edgeId(0, 8)] = 5.0;
    expected[G.edgeId(1, 3)] = 1.0;
    expected[G.edgeId(3, 7)] = 1.0;
    expected[G.edgeId(5, 6)] = 1.0;
    expected[G.edgeId(6, 7)] = 1.0;
    expected[G.edgeId(2, 4)] = 13.0;
    expected[G.edgeId(4, 9)] = 10.0;
    expected[G.edgeId(4, 8)] = 3.0;

    G.forEdges([&](node, node, edgeid eid) { EXPECT_DOUBLE_EQ(flow.get(eid), expected[eid]); });
    checkFlowConservation(G, flow);
}
} // namespace NetworKit
