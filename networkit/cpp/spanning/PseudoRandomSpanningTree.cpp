/*
 * PseudoRandomSpanningTree.cpp
 *
 *  Created on: 20.06.2015
 *      Author: Henning
 */

#include "PseudoRandomSpanningTree.h"
#include "../structures/UnionFind.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

PseudoRandomSpanningTree::PseudoRandomSpanningTree(const Graph& G): g(G) {

}

void PseudoRandomSpanningTree::run() {
	count n = g.numberOfNodes();
	Graph randTree(n);
	UnionFind uf(n);

	// sort edges in decreasing weight order
	std::vector<MyEdge> sortedEdges;
	g.forEdges([&](node u, node v, edgeweight ew) {
		double randVal = 1e-6 * (1.0 - 2.0 * Aux::Random::probability());

		MyEdge myEdge;
		myEdge.from = u;
		myEdge.to = v;
		myEdge.weight = ew + randVal;
		sortedEdges.push_back(myEdge);
	});
	std::sort(sortedEdges.begin(), sortedEdges.end());

	// process in decreasing weight order
	for (auto e: sortedEdges) {
		node u = e.from;
		node v = e.to;
//		INFO("process edge (", u, ", ", v, ") with weight ", e.weight);

		// if edge does not close cycle, add it to tree
		if (uf.find(u) != uf.find(v)) {
			randTree.addEdge(u, v);
			uf.merge(u, v);
		}
	}

	tree = randTree;
}

void PseudoRandomSpanningTree::runShuffle() {

	// TODO: handle disconnected graphs

	count n = g.numberOfNodes();
	Graph randTree(n);
	UnionFind uf(n);

	// prepare array to be shuffled
	std::vector<std::pair<node, node> > multEdges;
	g.forEdges([&](node u, node v) {
		count mult = 1; // g.degree(u) + g.degree(v);
		for (index i = 0; i < mult; ++i) {
			multEdges.push_back(std::make_pair(u, v));
		}
	});


	// shuffle edges to get their processing order
	std::random_shuffle(multEdges.begin(), multEdges.end());
	for (auto e: multEdges) {
		node u = e.first;
		node v = e.second;

		// if edge does not close cycle, add it to tree
		if (uf.find(u) != uf.find(v)) {
			randTree.addEdge(u, v);
			uf.merge(u, v);
		}
	}

	tree = randTree;
}

Graph PseudoRandomSpanningTree::getTree() {
	return tree;
}

} /* namespace NetworKit */
