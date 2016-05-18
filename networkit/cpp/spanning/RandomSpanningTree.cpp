/*
 * RandomSpanningTree.cpp
 *
 *  Created on: 20.06.2015
 *      Author: Henning
 */

#include "RandomSpanningTree.h"
#include "../structures/UnionFind.h"
#include "../graph/Sampling.h"
#include <random>

namespace NetworKit {

RandomSpanningTree::RandomSpanningTree(const Graph& G): g(G) {
	edges.clear();
	edges.resize(G.numberOfEdges());
	count i = 0;
	G.forEdges([&](node u, node v){
		edges[i] = std::make_pair(u,v);
		i ++;
	});
}

RandomSpanningTree::~RandomSpanningTree() {
}


void RandomSpanningTree::run() {

	// TODO: handle disconnected graphs

	count n = g.numberOfNodes();
	Graph randTree(n);
	count numVisited = 0;
	std::vector<bool> visited(n, false);

	// find and process root
	node curr = Sampling::randomNode(g);
	visited[curr] = true;
	numVisited++;

	while (numVisited < n) {
		// get random neighbor
		node neigh = g.randomNeighbor(curr);

		// if not seen before, insert tree edge
		if (! visited[neigh]) {
			randTree.addEdge(curr, neigh);
			visited[neigh] = true;
			++numVisited;
		}

		// move to neighbor
		curr = neigh;
	}

	tree = randTree;
}

void RandomSpanningTree::run2() {

	// TODO: handle disconnected graphs
	count n = g.numberOfNodes();
	// std::default_random_engine generator;
	// std::uniform_int_distribution<int> distribution(0,g.numberOfEdges()-1);
	std::random_shuffle (edges.begin(), edges.end());
	Graph randTree(n);
	UnionFind part(n);
	// INFO("Number of partitions: ", part.numberOfSubsets());
	count iter = 0, edgesTree = 0;
	while(edgesTree < n-1) {
	//	INFO("Number of partitions: ", part.numberOfSubsets());
		std::pair<node, node> rand_edge = edges[iter];
	//	INFO(rand_edge.first, " ", rand_edge.second);
		if (part.find(rand_edge.first) != part.find(rand_edge.second)) {
			// INFO("Partition of the first: ", rand_edge.first);
			// INFO("Partition of the second: ", rand_edge.second);
			// INFO("Merging two partitions. Before :", part.numberOfSubsets());
			randTree.addEdge(rand_edge.first, rand_edge.second);
			part.merge(part.find(rand_edge.first), part.find(rand_edge.second));
			edgesTree ++;
			// INFO("AFter: ", part.numberOfSubsets());
		}
		iter ++;
	}
	INFO("Iter: ", iter);
	assert(randTree.numberOfEdges() == n-1);
	tree = randTree;
}

Graph RandomSpanningTree::getTree() {
	return tree;
}

} /* namespace NetworKit */
