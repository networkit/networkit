/*
 * SpanningForest.cpp
 *
 *  Created on: 06.09.2015
 *      Author: Henning
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/SpanningForest.hpp>

namespace NetworKit {

SpanningForest::SpanningForest(const Graph& G): G(G) {

}

Graph SpanningForest::getForest() {
    return forest;
}

void SpanningForest::run() {
    forest = generate();
}

// CLS's code for generating spanning forest with BFS, please fixme!
Graph SpanningForest::generate() {
    Graph F = GraphTools::copyNodes(G);
    std::vector<bool> visited(G.upperNodeIdBound(), false);

    G.forNodes([&](node s){
        if (! visited[s]) {
            Traversal::BFSEdgesFrom(G, s, [&](node u, node v, edgeweight w, edgeid) {
                visited[u] = true;
                visited[v] = true;
                F.addEdge(u, v, w);
            });
        }
    });

    INFO("tree edges in SpanningForest: ", F.numberOfEdges());

    return F;

}

} /* namespace NetworKit */
