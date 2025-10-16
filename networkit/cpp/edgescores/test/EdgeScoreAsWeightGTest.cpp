/*
 * EdgeScoreGTest.cpp
 *
 *	Created on: 16.10.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 */
#include <gtest/gtest.h>

#include <networkit/edgescores/EdgeScoreAsWeight.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

static std::vector<double> makeSequentialScores(const Graph &G) {
    std::vector<double> s(G.numberOfEdges());
    for (edgeid eid = 0; eid < G.numberOfEdges(); ++eid) {
        s[eid] = static_cast<double>(eid + 1); // 1.0, 2.0, 3.0, ...
    }
    return s;
}

TEST(EdgeScoreAsWeight, testThrowsIfEdgesNotIndexed) {
    Graph G(2, /*weighted=*/false, /*directed=*/false);
    G.addEdge(0, 1);
    EdgeScoreAsWeight esw(G, std::vector<double>{}, /*squared=*/false, /*offset=*/0.0,
                          /*factor=*/1.0);
    EXPECT_THROW(
        {
            auto H = esw.calculate();
            (void)H;
        },
        std::runtime_error);
}

TEST(EdgeScoreAsWeight, testComputesWeightsLinear) {
    Graph G(3, /*weighted=*/false, /*directed=*/false);
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.indexEdges();

    const std::vector<double> scores = makeSequentialScores(G);

    constexpr double offset = 1.5;
    constexpr double factor = 2.0;

    EdgeScoreAsWeight esw(G, scores, /*squared=*/false, offset, factor);
    Graph H = esw.calculate();

    EXPECT_EQ(H.numberOfNodes(), G.numberOfNodes());
    EXPECT_EQ(H.numberOfEdges(), G.numberOfEdges());
    EXPECT_TRUE(H.isWeighted());

    G.forEdges([&](node u, node v) {
        EXPECT_DOUBLE_EQ(H.weight(u, v), offset + factor * scores[G.edgeId(u, v)]);
    });
}

TEST(EdgeScoreAsWeight, testComputesWeightsSquared) {
    Graph G(4, /*weighted=*/false, /*directed=*/false, /*edgesIndex*/ true);
    G.addEdge(0, 1);
    G.addEdge(0, 2);
    G.addEdge(1, 3);

    const std::vector<double> scores = makeSequentialScores(G);
    constexpr double offset = -3.0;
    constexpr double factor = 0.5;

    EdgeScoreAsWeight esw(G, scores, /*squared=*/true, offset, factor);
    Graph H = esw.calculate();

    EXPECT_EQ(H.numberOfEdges(), G.numberOfEdges());
    EXPECT_TRUE(H.isWeighted());

    G.forEdges([&](node u, node v) {
        EXPECT_DOUBLE_EQ(H.weight(u, v),
                         offset + factor * scores[G.edgeId(u, v)] * scores[G.edgeId(u, v)]);
    });
}

} // namespace NetworKit
