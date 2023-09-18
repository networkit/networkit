/*
* Bipartit.cpp
*
* Created on: 18.09.2023
*     Author: Michael Kaibel
 */

#include <gtest/gtest.h>

#include "networkit/bipartit/Bipartit.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/generators/ErdosRenyiGenerator.hpp"

namespace NetworKit {

class BipartitGTest : public testing::Test {};

bool bipartitPartition(const Graph &G, const Partition &partition) {
    for (node v : G.nodeRange()) {
        if (partition[v] == none)
            return false;

        assert(partition[v] == 0 or partition[v] == 1);

        for (node w : G.neighborRange(v)) {
            if (v <= w)
                continue;

            if (partition[v] == partition[w])
                return false;
        }
    }

    return true;
}

TEST_F(BipartitGTest, testBipartitTinyBipartit) {
    Graph G(8);
    G.addEdge(0, 5);
    G.addEdge(0, 6);
    G.addEdge(1, 6);
    G.addEdge(1, 7);
    G.addEdge(2, 4);
    G.addEdge(2, 5);
    G.addEdge(2, 6);
    G.addEdge(2, 7);
    G.addEdge(3, 7);

    Bipartit bipartit(G);
    bipartit.run();

    EXPECT_TRUE(bipartit.isBipartit());

    EXPECT_TRUE(bipartitPartition(G, bipartit.getPartition()));
}

TEST_F(BipartitGTest, testBipartitTinyOddCircle) {
    Graph G(8);
    G.addEdge(0, 2);
    G.addEdge(2, 5);
    G.addEdge(5, 4);
    G.addEdge(4, 1);
    G.addEdge(1, 0);
    G.addEdge(5, 6);
    G.addEdge(6, 7);
    G.addEdge(7, 0);

    Bipartit bipartit(G);
    bipartit.run();

    EXPECT_FALSE(bipartit.isBipartit());

    EXPECT_ANY_THROW(bipartit.getPartition());
}

TEST_F(BipartitGTest, testBipartitLargeRandom) {
    Aux::Random::setSeed(42, false);
    Graph G = ErdosRenyiGenerator(200, 0.01, false).generate();

    Bipartit bipartit(G);
    bipartit.run();

    if (bipartit.isBipartit()) {
        EXPECT_TRUE(bipartitPartition(G, bipartit.getPartition()));
    }
}

} // NetworKit