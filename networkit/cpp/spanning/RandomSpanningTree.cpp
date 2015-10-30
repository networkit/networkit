/*
 * RandomSpanningTree.cpp
 *
 *  Created on: 20.06.2015
 *      Author: Henning
 */

#include "RandomSpanningTree.h"
#include "../graph/Sampling.h"

namespace NetworKit {

RandomSpanningTree::RandomSpanningTree(const Graph& G): g(G) {
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

Graph RandomSpanningTree::getTree() {
	return tree;
}

} /* namespace NetworKit */
