/*
* WattsStrogatzGenerator.cpp
*
*  Created on: 09.07.2014
*      Author: Simon Bischof
*/

// networkit-format

#include <limits>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/generators/RegularRingLatticeGenerator.hpp>
#include <networkit/generators/WattsStrogatzGenerator.hpp>

namespace NetworKit {

WattsStrogatzGenerator::WattsStrogatzGenerator(count nNodes, count nNeighbors, double p)
    : nNodes(nNodes), nNeighbors(nNeighbors), p(p) {
		if (nNeighbors*2 >= nNodes-1) {
			throw std::runtime_error("nNeighbors*2 cannot be equal to nNodes-1.");
		} else if (nNeighbors >= nNodes / 2 - 1) {
            nNeighbors = nNodes / 2 - 1;
        }
	}

Graph WattsStrogatzGenerator::generate() {
    //generate regular ring lattice as initial graph
    RegularRingLatticeGenerator R(nNodes, nNeighbors);
    Graph G = R.generate();
    std::uniform_int_distribution<node> dist(0, nNodes - 1);
    std::mt19937_64 generator(Aux::Random::integer(0, std::numeric_limits<count>::max()));
    std::uniform_real_distribution<double> prob{0.0, std::nexttoward(1.0, 2.0)};

    //list all valid end nodes for rewiring an edge
    auto validEndNode = [&G](node u, node v) {
        return u != v && !G.hasEdge(u, v);
    };

    std::vector<node> prevNeighbors;
    prevNeighbors.reserve(GraphTools::maxDegree(G));

    // rewire according to Watts and Strogatz model
    for (const auto u : G.nodeRange()) {
        if (G.degree(u) == (nNodes - 1))
            continue;
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
            if (prob(generator) < p) {
                bool found = false;
                node w = none;
                while(!found) {
                    w = dist(generator);
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
