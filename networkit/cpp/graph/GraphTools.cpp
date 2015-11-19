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

std::vector<node> invertContinuousNodeIds(std::unordered_map<node,node>& nodeIdMap, const Graph& G) {
	assert(nodeIdMap.size() == G.numberOfNodes());
	std::vector<node> invertedIdMap(G.numberOfNodes() + 1);
	// store upper node id bound
	invertedIdMap[G.numberOfNodes()] = G.upperNodeIdBound();
	// inverted node mapping
	for (auto& x : nodeIdMap) {
		invertedIdMap[x.second] = x.first;
	}
	return invertedIdMap;
}

Graph restoreGraph(std::vector<node>& invertedIdMap, const Graph& G) {
	// with the inverted id map and the compacted graph, generate the original graph again
	Graph Goriginal(invertedIdMap[invertedIdMap.size()-1],G.isWeighted(), G.isDirected());
	index current = 0;
	Goriginal.forNodes([&](node u){
		if (invertedIdMap[current] == u) {
			G.forNeighborsOf(current,[&](node v){
				Goriginal.addEdge(u,invertedIdMap[v]);
			});
			++current;
		} else {
			Goriginal.removeNode(u);
		}
	});
	return Goriginal;
}

}

}