#include <algorithm>
#include <random>
#include <unordered_map>

#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

namespace GraphTools {

Graph getCompactedGraph(const Graph& graph, const std::unordered_map<node,node>& nodeIdMap) {
	return getRemappedGraph(graph, nodeIdMap.size(), [&] (node u) {
	    const auto it = nodeIdMap.find(u);
	    assert(it != nodeIdMap.cend());
	    return it->second;
	});
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

std::vector<node> invertContinuousNodeIds(const std::unordered_map<node,node>& nodeIdMap, const Graph& G) {
	assert(nodeIdMap.size() == G.numberOfNodes());
	std::vector<node> invertedIdMap(G.numberOfNodes() + 1);
	// store upper node id bound
	invertedIdMap[G.numberOfNodes()] = G.upperNodeIdBound();
	// inverted node mapping
	for (const auto x : nodeIdMap) {
		invertedIdMap[x.second] = x.first;
	}
	return invertedIdMap;
}

Graph restoreGraph(const std::vector<node>& invertedIdMap, const Graph& G) {
	// with the inverted id map and the compacted graph, generate the original graph again
	Graph Goriginal(invertedIdMap.back(), G.isWeighted(), G.isDirected());
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

Graph extractLargestConnectedComponent(const Graph &G, bool compactGraph) {
	if (!G.numberOfNodes())
		return G;

	ConnectedComponents cc(G);
	cc.run();

	auto compSizes = cc.getComponentSizes();
	if (compSizes.size() == 1)
			return G;

	auto largestCC = std::max_element(
		compSizes.begin(), compSizes.end(),
		[](const std::pair<index, count> &x, const std::pair<index, count> &y) {
			return x.second < y.second;
		});

	std::unordered_map<node, node> continuousNodeIds;
	index nextId = 0;
	G.forNodes([&](node u) {
		if (cc.componentOfNode(u) == largestCC->first)
			continuousNodeIds[u] = nextId++;
	});

	if (compactGraph) {
		Graph S = getRemappedGraph(G, largestCC->second,
			[&](node u) { return continuousNodeIds[u]; },
			[&](node u) { return cc.componentOfNode(u) != largestCC->first; });
		return S;
	} else {
		Graph S(G);
		auto components = cc.getComponents();
		for (count i = 0; i < components.size(); ++i) {
			if (i != largestCC->first) {
				for (node u : components[i])
					S.removeNode(u);
			}
		}
		return S;
	}
}

} // namespace GraphTools
} // namespace NetworKit
