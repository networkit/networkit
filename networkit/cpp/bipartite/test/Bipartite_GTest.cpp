/*
 * BipartitGTest.cpp
 *
 * Created on: 18.09.2023
 *     Author: Michael Kaibel
 */

#include <gtest/gtest.h>

#include <random>

#include "networkit/bipartite/Bipartite.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/generators/ErdosRenyiGenerator.hpp"

namespace NetworKit {

class BipartitGTest : public testing::Test {};

bool isBipartitePartition(const Graph &G, const Partition &partition) {
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

bool isOddCircle(const Graph &G, const std::vector<node> &path) {
    if (path.size() % 2 == 0)
        return false;
    for (index i = 0; i < path.size(); i++)
        for (index j = i+1; j < path.size(); j++)
            if (path[i] == path[j])
                throw std::runtime_error("Not a path");

    for (index i = 0; i < path.size(); i++)
        if (not G.hasEdge(path[i], path[(i + 1) % path.size()]))
            return false;

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

    Bipartite bipartite(G);
    bipartite.run();

    EXPECT_TRUE(bipartite.isBipartite());

    EXPECT_TRUE(isBipartitePartition(G, bipartite.getPartition()));

    EXPECT_ANY_THROW(bipartite.getOddCycle());
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

    Bipartite bipartite(G);
    bipartite.run();

    EXPECT_FALSE(bipartite.isBipartite());

    EXPECT_ANY_THROW(bipartite.getPartition());

    EXPECT_TRUE(isOddCircle(G, bipartite.getOddCycle()));
}

TEST_F(BipartitGTest, testBipartitLargeRandom) {
    Aux::Random::setSeed(42, false);
    Graph G = ErdosRenyiGenerator(200, 0.01, false).generate();

    Bipartite bipartite(G);
    bipartite.run();

    if (bipartite.isBipartite()) {
        EXPECT_TRUE(isBipartitePartition(G, bipartite.getPartition()));

        EXPECT_ANY_THROW(bipartite.getOddCycle());
    } else {
        EXPECT_ANY_THROW(bipartite.getPartition());

        EXPECT_TRUE(isOddCircle(G, bipartite.getOddCycle()));
    }
}

//Since a random ErdosRenyi-Graph almost certainly has an odd circle we also do a test with an explicitly bipartit graph
TEST_F(BipartitGTest, testBipartitLargeRandomBipartit) {
    node n = 200;

    std::mt19937 rng(42);
    std::uniform_int_distribution<node> dis(1000000);

    std::vector<node> nodes(n);
    for (index i = 0; i < n; i++)
        nodes[i] = i;

    for (index i = 0; i < n; i++)
        std::swap(nodes[i], nodes[i + (dis(rng) % (n - i))]);

    Graph G(n);

    for (index i = 0; i < 5*n; i++) {
        node v = nodes[dis(rng) % (n/2)], w = nodes[(n/2) + (dis(rng) % (n/2))];

        if (not G.hasEdge(v, w))
            G.addEdge(v, w);
    }

    Bipartite bipartite(G);
    bipartite.run();

    EXPECT_TRUE(bipartite.isBipartite());

    EXPECT_TRUE(isBipartitePartition(G, bipartite.getPartition()));

    EXPECT_ANY_THROW(bipartite.getOddCycle());
}

} // namespace NetworKit
