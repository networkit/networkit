/*
 * TriangleScoreGTest.cpp
 *
 *  Created on: 05.09.2023
 *      Author: Angus Petholz
 */

#include <gtest/gtest.h>

#include <networkit/sparsification/ChanceCorrectedTriangleScore.hpp>
#include <networkit/sparsification/SCANStructuralSimilarityScore.hpp>

namespace NetworKit {

class TriangleScoreGTest : public testing::Test {};

Graph initGraph() {
    Graph G(5);
    G.addEdge(0, 1); // 0     G:  0 - 1
    G.addEdge(0, 2); // 1         | / |
    G.addEdge(1, 2); // 2         2 - 3
    G.addEdge(1, 3); // 3         |
    G.addEdge(2, 3); // 4         4
    G.addEdge(2, 4); // 5
    G.indexEdges();  // EdgeId
    return G;
}

TEST_F(TriangleScoreGTest, testChanceCorrectedTriangleScore) {
    Graph G = initGraph();

    // edges 0,1,3,4 are contained in one triangle, edge 2 is contained in two triangles
    std::vector<count> triangles = {1, 1, 2, 1, 1, 0};
    std::vector<double> expectedRes = {1.5, 1, 1, 1.5, 1, 1};

    ChanceCorrectedTriangleScore cct(G, triangles);
    cct.run();
    auto res = cct.scores();

    EXPECT_EQ(res.size(), G.numberOfEdges());
    for (index i = 0; i < res.size(); ++i)
        EXPECT_EQ(res[i], expectedRes[i]);
}

TEST_F(TriangleScoreGTest, testSCANStructuralSimiliarityScore) {
    Graph G = initGraph();

    std::vector<count> triangles = {1, 1, 2, 1, 1, 0};
    std::vector<double> expectedRes = {0.57, 0.51, 0.67, 0.57, 0.51, 0.31};

    SCANStructuralSimilarityScore ssss(G, triangles);
    ssss.run();
    auto res = ssss.scores();

    EXPECT_EQ(res.size(), G.numberOfEdges());
    for (index i = 0; i < res.size(); ++i)
        EXPECT_NEAR(res[i], expectedRes[i], 0.01);
}

} // namespace NetworKit
/* namespace NetworKit */
