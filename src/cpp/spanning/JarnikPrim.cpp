/*
 * JarnikPrim.cpp
 *
 *  Created on: 13.05.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "JarnikPrim.h"

namespace NetworKit {

std::vector<std::vector<Edge>> JarnikPrim::run(const Graph &G) {
	std::vector<std::vector<Edge>> msf;
	std::vector<Edge> mst;
	std::vector<bool> visited(G.numberOfNodes(), false);
	node root = 0;
	visited[root] = true;
	count visitedNodes = 1; // root already visited

	// build priority queue and add adjacent edges
	std::map<index, Edge> indexToEdge;

	index edgeIdx = 0;
	std::vector<std::pair<double, index>> rootNeighbors;
	// add edges of root to rootNeighbors
	G.forWeightedNeighborsOf(root, [&](node v, double weight) {
		Edge edge = std::make_pair(root, v);
		rootNeighbors.push_back(std::make_pair(weight, edgeIdx));
		indexToEdge.insert(std::make_pair(edgeIdx++, edge));
	});

	Aux::PrioQueue<double, index> pq(rootNeighbors);
	while (visitedNodes < visited.size()) {
		if (pq.size() > 0) {
			Edge edge = indexToEdge[pq.extractMin().second];
			if (!visited[edge.second]) { // not visited neighbor
				visited[edge.second] = true;
				visitedNodes++;
				mst.push_back(edge);

				G.forWeightedNeighborsOf(edge.second, [&](node v, double weight){
					if (v != edge.first) { // do not insert the current edge
						pq.insert(weight, edgeIdx);
						indexToEdge.insert(std::make_pair(edgeIdx++, std::make_pair(edge.second, v)));
					}
				});
			}
		} else { // create new mst and save the current one
			msf.push_back(mst);
			mst = std::vector<Edge>();

			// add new root node for this connected component
			bool foundUnvisitedNode = false;
			G.forNodes([&](){return !foundUnvisitedNode;}, [&](node u){
				if (!visited[u]) {
					visited[u] = true;
					visitedNodes++;
					foundUnvisitedNode = true;

					// add neighbors to priority queue
					G.forWeightedNeighborsOf(u, [&](node v, double weight) {
						Edge edge = std::make_pair(u, v);
						pq.insert(weight, edgeIdx);
						indexToEdge.insert(std::make_pair(edgeIdx++, edge));
					});
				}
			});
		}
	}

	msf.push_back(mst);

	return msf;
}

} /* namespace NetworKit */
