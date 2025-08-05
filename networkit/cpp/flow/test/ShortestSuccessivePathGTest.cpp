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

TEST_F(ShortestSuccessivePathGTest, constructorSucceedsWhenValid) {
    Graph G(3, /*weighted*/ true, /*directed*/ true);
    G.indexEdges();
    auto eAttr = G.attachEdgeDoubleAttribute(capacityName);
    auto nAttr = G.attachNodeDoubleAttribute(supplyName);

    // populate attributes
    G.forNodes([&](node u) { nAttr.set(u, 0.0); });
    G.addEdge(0, 1, /*weight=*/1.5, /*checkMulti=*/false);
    G.addEdge(1, 2, /*weight=*/2.5, /*checkMulti=*/false);

    eAttr.set(G.edgeId(0, 1), 10.0);
    eAttr.set(G.edgeId(1, 2), 20.0);

    EXPECT_NO_THROW({ MinFlowShortestSuccessivePath test(G, capacityName, supplyName); });
}

} // namespace NetworKit
