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
    // edges are not indexed by default
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
    // attach only the node attribute
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
    // attach only the edge attribute
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
    // Build a minimal valid directed, weighted, indexed graph
    Graph G(2, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();

    // Attach the capacity and supply attributes
    auto capacities   = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    // Add one edge and set a negative capacity on it
    G.addEdge(0, 1, /*weight=*/1.0);
    auto eid = G.edgeId(0, 1);
    capacities.set(eid, -5.0);  // invalid negative capacity

    // Supplies can all be zero
    G.forNodes([&](node u) {
        supply.set(u, 0.0);
    });

    // Expect the constructor to catch the negative capacity
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

    // populate attributes
    G.forNodes([&](node u) { supply.set(u, 0.0); });
    G.addEdge(0, 1, /*cost=*/1.5);
    G.addEdge(1, 2, /*cost=*/2.5);

    capacities.set(G.edgeId(0, 1), 10.0);
    capacities.set(G.edgeId(1, 2), 20.0);

    EXPECT_NO_THROW({ MinFlowShortestSuccessivePath test(G, capacityName, supplyName); });
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

    // 7) Check that total cost = 3*2 + 2*3 = 12
    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 12.0);
}

TEST_F(ShortestSuccessivePathGTest, ComplexMultiSourceNegativeCost) {
    // 1) Build a directed, weighted graph on 6 nodes
    Graph G(6, /*weighted=*/true, /*directed=*/true);
    G.indexEdges();

    // 2) Attach capacity & supply attributes
    auto capacities = G.attachEdgeDoubleAttribute(capacityName);
    auto supply = G.attachNodeDoubleAttribute(supplyName);

    // 3) Set up supplies: node 0 → +4, node 1 → +3, node 5 → -7
    for (node u = 0; u < 6; ++u) {
        supply.set(u, 0.0);
    }
    supply.set(0, +4.0);
    supply.set(1, +3.0);
    supply.set(5, -7.0);

    auto addEdgeAndCapacity = [&](node u, node v, double cost, double capacity) {
        G.addEdge(u, v, cost);
        capacities.set(G.edgeId(u, v), capacity);
    };

    addEdgeAndCapacity(0, 2, /*cost=*/2.0, /*capacity=*/4.0);
    addEdgeAndCapacity(2, 5, /*cost=*/1.0, /*capacity=*/4.0);
    addEdgeAndCapacity(0, 3, /*cost=*/1.0, /*capacity=*/2.0);
    addEdgeAndCapacity(3, 5, /*cost=*/3.0, /*capacity=*/2.0);
    addEdgeAndCapacity(1, 3, /*cost=*/-1.0,/*capacity=*/3.0);
    addEdgeAndCapacity(3, 4, /*cost=*/1.0, /*capacity=*/3.0);
    addEdgeAndCapacity(4, 5, /*cost=*/2.0, /*capacity=*/3.0);
    addEdgeAndCapacity(1, 2, /*cost=*/3.0, /*capacity=*/3.0);

    MinFlowShortestSuccessivePath solver(G, capacityName, supplyName);
    solver.run();

    EXPECT_DOUBLE_EQ(solver.getTotalCost(), 18.0);
}

} // namespace NetworKit
