/*
 * CoreDecomposition.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Henning
 */

#include "CoreDecomposition.h"

namespace NetworKit {

CoreDecomposition::CoreDecomposition() {

}

CoreDecomposition::~CoreDecomposition() {

}

std::vector<count> CoreDecomposition::run(const Graph& G) {
#if 0 // list does not compile after pull and merge with Brueckner code

    using namespace std;

	vector<count> degree(G.numberOfNodes());
    vector<list<node>> nodeLists(G.numberOfNodes(), list<node>());
    vector<pair<list<node>, list<node>::iterator>> nodeLocation(G.numberOfNodes());

    G.forNodes([&](node v) {
        degree[v] = G.degree(v);
        nodeLists[degree[v]].push_front(v);
        nodeLocation[v] = pair<list<node>, list<node>::iterator>(nodeLists[degree[v]], nodeLists[degree[v]].begin());
    });

    vector<count> coreness(G.numberOfNodes());
    Graph G2 = G;
    index i = 0;

    while (G2.numberOfNodes() > 0) {
        for (node u : nodeLists[i]) {
            coreness[u] = i;
            G2.forNeighborsOf(u, [&](node v) {
                if (--degree[v] >= i) {
                    nodeLocation[v].first.erase(nodeLocation[v].second);
                    nodeLists[degree[v]].push_back(v);
                    nodeLocation[v] = pair<list<node>, list<node>::iterator>(nodeLists[degree[v]], prev(nodeLists[degree[v]].end()));
                }
                G2.removeEdge(u, v);
            });
            G2.removeNode(u);
        }
        i++;
    }

    return coreness;
#endif
}

} /* namespace NetworKit */
