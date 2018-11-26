/*
 * GlobalCurveballGTest.cpp
 *
 *  Created on: 24.05.2018
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#include <gtest/gtest.h>

#include "../GlobalCurveball.h"
#include "../../graph/Graph.h"
#include "../../generators/ErdosRenyiGenerator.h"
#include "../../generators/HyperbolicGenerator.h"

namespace NetworKit {

class GlobalCurveballGTest : public ::testing::Test  {
protected:
    void checkWithGraph(Graph&);
};

void GlobalCurveballGTest::checkWithGraph(Graph& G) {
    node numNodes = G.numberOfNodes();
    const count numTrades = 5;

    std::vector<node> degrees(numNodes + 1);

    // Add edge to node 0, if isolated node
    // If 0 itself is isolated, add new node and connect 0 to it
    G.forNodes([&](node u) {
        if (G.degree(u) > 0)
            degrees[u] = G.degree(u);
        else {
            if (u == 0) {
                numNodes++;
                G.addEdge(0, numNodes - 1);
                degrees[0]++;
                degrees[numNodes - 1] = 1;
            } else {
                G.addEdge(u, 0);
                degrees[0]++;
                degrees[u] = 1;
            }
        }
    });


    GlobalCurveball algo(G, numTrades);
    algo.run();

    // check degrees
    Graph outG = algo.getGraph();
    outG.forNodes([&](node u){
        ASSERT_EQ(degrees[u], outG.degree(u));
    });
}



TEST_F(GlobalCurveballGTest, testCurveballErdosRenyi) {
    Aux::Random::setSeed(1, false);

    node numNodes = 1000;
    ErdosRenyiGenerator generator(numNodes, 0.01);
    Graph G = generator.generate();

    this->checkWithGraph(G);
}

TEST_F(GlobalCurveballGTest, testCurveballHyperbolic) {
    Aux::Random::setSeed(1, false);

    node numNodes = 1000;
    HyperbolicGenerator generator(numNodes);
    Graph G = generator.generate();

    this->checkWithGraph(G);
}

} // namespace NetworKit
