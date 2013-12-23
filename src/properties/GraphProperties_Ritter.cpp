/*
 * GraphProperties_Ritter.cpp
 *
 *  Created on: 24.11.2013
 *      Author: mrittre
 */

#include "GraphProperties_Ritter.h"
#include "../graph/BFS.h"

namespace NetworKit {

const count GraphProperties_Ritter::INF_DIST = std::numeric_limits<count>::max();

GraphProperties_Ritter::GraphProperties_Ritter() {

}

GraphProperties_Ritter::~GraphProperties_Ritter() {

}

count GraphProperties_Ritter::diameter(const Graph& G) {
	count max_ecc = 0;

	// the diameter is the maximum eccentricity of all nodes
	G.forNodes([&] (node v) {
		count ecc = eccentricity(G, v);
		if (ecc > max_ecc) {
			max_ecc = ecc;
		}
	});

	return max_ecc;
}

count GraphProperties_Ritter::diameterOfTree(const Graph& T, node root) {
	std::pair<count, node> ecc = eccentricityWithNode(T, root);
	return eccentricity(T, ecc.second);
}

count GraphProperties_Ritter::eccentricity(const Graph& G, node v) {
	return eccentricityWithNode(G, v).first;
}

std::pair<count, index> GraphProperties_Ritter::eccentricityWithNode(const Graph& G, node v) {
	const count n = G.numberOfNodes();
	/*
	 * 2 possible ways:
	 * a) use BFS class and search for index with max distance => runtime = 2 n
	 * b) do BFS our self and remember last index (also one with highest distance) => runtime = n
 	 */
	// a)
	// std::vector<count> distances = BFS().run(G, v);
	// index maxDistanceIndex = 0;
	// for (index i = 1; i < n; i++) {
	// 	if (distances[i] > distances[maxDistanceIndex]) {
	// 		maxDistanceIndex = i;
	// 	}
	// }
	// return std::make_pair(distances[maxDistanceIndex], maxDistanceIndex);

	// b)
	std::vector<count> distances(n, INF_DIST);
	std::queue<node> que;

	distances[v] = 0;
	que.push(v);

	node lastNode;
	count processedNodes = 0;
	while (!que.empty()) {
		processedNodes++;
		node current = que.front();
		que.pop();

		// insert untouched neighbors into queue
		G.forNeighborsOf(current, [&](node neighbor) {
			if (distances[neighbor] == INF_DIST) {
				que.push(neighbor);
				distances[neighbor] = distances[current] + 1;
			}
		});

		lastNode = current;
	}

	if (processedNodes == n) {
		// all nodes were visited => G is connected
		return std::make_pair(distances[lastNode], lastNode);
	} else {
		// at least one node was not visited and has distance = infinity (G is not connected)
		return std::make_pair(INF_DIST, 0);
	}
}

Graph GraphProperties_Ritter::spanningTree(const Graph& G, node root) {
	const count n = G.numberOfNodes();
	Graph T(n);
	
	std::vector<bool> inSpanningTree(n, false);
	inSpanningTree[root] = true;

	std::queue<node> q;
	q.push(root);

	while (!q.empty()) {
		node current = q.front();
		q.pop();

		G.forNeighborsOf(current, [&](node v) {
			if (!inSpanningTree[v]) {
				q.push(v);
				T.addEdge(current, v);
				inSpanningTree[v] = true;
			}
		});
	}

	return T;
}

} /* namespace NetworKit */
