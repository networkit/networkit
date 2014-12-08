#include "GraphTools.h"
#include <unordered_map>
#include "../graph/Graph.h"

namespace NetworKit {

namespace GraphTools {

Graph getCompactedGraph(const Graph& graph) {
	auto nodeIdMap = getContinuousNodeIds(graph);
	Graph Gcompact(nodeIdMap.size(),graph.isWeighted(),graph.isDirected());
	auto copyEdge = [&Gcompact,&nodeIdMap](node u, node v, edgeweight ew) {
		Gcompact.addEdge(nodeIdMap[u], nodeIdMap[v], ew);
	};
	graph.forEdges(copyEdge);
	return Gcompact;
}

std::unordered_map<node,node> getContinuousNodeIds(const Graph& graph) {
	std::unordered_map<node,node> nodeIdMap;
	count continuousId = 0;
	auto addToMap = [&nodeIdMap,&continuousId](node v) {
		nodeIdMap.insert(std::make_pair(v,continuousId++));
	};
	graph.forNodes(addToMap);
	return nodeIdMap;
}

Graph toUndirected(const Graph& graph) {
	return Graph(graph,graph.isWeighted(),false);
}

}

}