/*
* BfsSpanningForest.cpp
*
*  Created on: Aug 7, 2014
*      Author: Christian Staudt
*/
#include "BfsSpanningForest.h"

namespace NetworKit {


BfsSpanningForest::BfsSpanningForest(const Graph& G) : G(G) {

}

Graph BfsSpanningForest::generate() {
	Graph F = G.copyNodes();
	std::vector<bool> visited(G.upperNodeIdBound(), false);

	G.forNodes([&](node s){
		if (! visited[s]) {
			G.BFSEdgesFrom(s, [&](node u, node v, edgeweight w, edgeid eid) {
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
