/*
 * AssortativityGTest.cpp
 *
 *  Created on: 04.03.2022
 *      Author: Lucas Petersen
 */

#include <gtest/gtest.h>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/correlation/Assortativity.hpp>
#include <networkit/io/METISGraphReader.hpp>

namespace NetworKit {

class AssortativityGTest : public testing::Test {};

TEST_F(AssortativityGTest, testAssortivityPartitionSmallGraph) {
    /* Graph:
       0   3   7
        \ /|\ /
         2 5 6
        / \|/ \
       1   4   8
    */
    Graph G(9);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(3, 6);
    G.addEdge(4, 5);
    G.addEdge(4, 6);
    G.addEdge(6, 7);
    G.addEdge(6, 8);

    /*
    Count both orientations of every undirected edge. The mixing matrix is
    normalized by 2m:
    |0.2 0.1 0.0|
    |0.1 0.2 0.1|
    |0.0 0.1 0.2|
    Its trace is 0.6 and the product of its marginals sums to 0.34.
    */
    std::vector<index> data = {0, 0, 0, 1, 1, 1, 2, 2, 2};
    Partition P(data);
    double coeff = 13. / 33.;

    Assortativity assort(G, P);
    assort.run();
    EXPECT_NEAR(coeff, assort.getCoefficient(), 1e-8);
}

TEST_F(AssortativityGTest, testAssortivityAttributeSmallGraph) {
    /* Graph:
       0   3   6
        \ / \ /
         2   5
        / \ / \
       1   4   7
    */
    Graph G(8);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.addEdge(5, 6);
    G.addEdge(5, 7);

    /*
    In this case it is equivalent to the Pearson correlation between:
    X = [1.0, 1.1, 1.2, 1.2, 2.0, 2.1, 3.0, 3.0]
    and
    Y = [1.2, 1.2, 2.0, 2.1, 3.0, 3.0, 3.1, 3.2]
    */
    std::vector<double> attribute = {1.0, 1.1, 1.2, 2.0, 2.1, 3.0, 3.1, 3.2};
    double coeff = 0.88237347;

    Assortativity assort(G, attribute);
    assort.run();
    EXPECT_NEAR(coeff, assort.getCoefficient(), 1e-8);
}

TEST_F(AssortativityGTest, testNominalAssortativityEdgeOrientations) {
    for (bool directed : {false, true}) {
        Graph G(4, false, directed);
        G.addEdge(0, 1);
        G.addEdge(1, 2);
        G.addEdge(2, 3);
        Partition groups(std::vector<index>{0, 0, 1, 1});
        Assortativity assort(G, groups);
        assort.run();
        // Undirected: [[1/3, 1/6], [1/6, 1/3]], giving 1/3.
        // Directed: [[1/3, 1/3], [0, 1/3]], giving 2/5.
        EXPECT_NEAR(assort.getCoefficient(), directed ? 2. / 5. : 1. / 3., 1e-12);
    }
}

TEST_F(AssortativityGTest, testNominalAssortativitySelfLoops) {
    for (bool directed : {false, true}) {
        Graph G(2, false, directed);
        G.addEdge(0, 0);
        G.addEdge(0, 1);
        G.addEdge(1, 1);
        Partition groups(std::vector<index>{0, 1});
        Assortativity assort(G, groups);
        assort.run();
        // A self-loop contributes once; other undirected edges contribute
        // once in each direction. Thus every undirected matrix entry is 1/4.
        EXPECT_NEAR(assort.getCoefficient(), directed ? 2. / 5. : 0., 1e-12);
    }
}

} // namespace NetworKit
