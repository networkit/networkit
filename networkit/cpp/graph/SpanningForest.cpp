/*
* SpanningForest.cpp
*
*  Created on: Aug 7, 2014
*      Author: Christian Staudt
*/
#include "SpanningForest.h"

namespace NetworKit {


SpanningForest::SpanningForest(const Graph& G) : G(G) {

}

Graph SpanningForest::generate() {
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

	return F;
}

} /* namespace NetworKit */
