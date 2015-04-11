#include "GraphTools.h"
#include <unordered_map>
#include "../graph/Graph.h"
#include <random>

namespace NetworKit {

namespace GraphTools {

Graph getCompactedGraph(const Graph& graph, std::unordered_map<node,node>& nodeIdMap) {
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

std::unordered_map<node,node> getRandomContinuousNodeIds(const Graph& graph) {
	std::unordered_map<node,node> nodeIdMap;
	std::vector<node> nodes;
	nodes.reserve(graph.numberOfNodes());

	graph.forNodes([&](node u) {
		nodes.push_back(u);
	});

	std::shuffle(nodes.begin(), nodes.end(), Aux::Random::getURNG());

	count continuousId = 0;
	for (node v : nodes) {
		nodeIdMap.insert(std::make_pair(v,continuousId++));
	};

	return nodeIdMap;
}

}

}