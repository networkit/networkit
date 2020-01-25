/*
* WattsStrogatzGenerator.cpp
*
*  Created on: 09.07.2014
*      Author: Simon Bischof
*/

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/RegularRingLatticeGenerator.hpp>
#include <networkit/generators/WattsStrogatzGenerator.hpp>

namespace NetworKit {

WattsStrogatzGenerator::WattsStrogatzGenerator(count nNodes, count nNeighbors, double p)
    : nNodes(nNodes), nNeighbors(nNeighbors), p(p) {}

Graph WattsStrogatzGenerator::generate() {
    //generate regular ring lattice as initial graph
    RegularRingLatticeGenerator R(nNodes, nNeighbors);
    Graph G = R.generate();

    //list all valid end nodes for rewiring an edge
    auto validEndNode = [&](node u, node v) {
        return ! G.hasEdge(u, v) && u != v;
    };

    //rewire according to Watts and Strogatz model
    for (node u = 0; u < nNodes; u++) {
        /* save edges before rewiring incident edges for a node
         * because they may get changed through rewiring
         */
        std::vector<node> prevNeighbors;
        G.forNeighborsOf(u, [&](node v) {
            if (u < v) {
                prevNeighbors.push_back(v);
            }
        });
        for (node v : prevNeighbors) {
            if (Aux::Random::probability() < p) {
                bool found = false;
                node w;
                while(!found) {
                    w = Aux::Random::integer(nNodes-1);
                    if (validEndNode(u,w)) {
                        found = true;
                    }
                }
                G.removeEdge(u, v);
                G.addEdge(u, w);
            }
        }
    }

    return G;
}

} /* namespace NetworKit */
