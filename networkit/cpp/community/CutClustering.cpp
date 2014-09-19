/*
 * Author: Michael Hamann <michael.hamann@kit.edu>
 */

#include "CutClustering.h"
#include "../flow/EdmondsKarp.h"

#include <sstream>

NetworKit::CutClustering::CutClustering(NetworKit::edgeweight alpha) : alpha(alpha) { }

NetworKit::Partition NetworKit::CutClustering::run(const NetworKit::Graph &G) {
	Partition result(G.upperNodeIdBound());
	result.setUpperBound(G.upperNodeIdBound());
	
	Graph graph(G, true, false);
	
	node t = graph.addNode();
	
	graph.forNodes([&](node u) {
		if (u != t) {
			graph.addEdge(u, t, alpha);
		}
	});
	
	graph.indexEdges();
	
	EdmondsKarp flowAlgo;
	
	graph.forNodes([&](node u) {
		if (u != t && !result.contains(u)) {
			std::vector<node> sourceSet;
			flowAlgo.run(graph, u, t, sourceSet);
			
			for (node v : sourceSet) {
				result[v] = u;
			}
		}
	});
	
	return result;
}

std::string NetworKit::CutClustering::toString() const {
	std::stringstream stream;
	
	stream << "CutClustering(" << alpha << ")";
	return stream.str();
}

