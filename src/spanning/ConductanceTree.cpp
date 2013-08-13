/*
 * ConductanceTree.cpp
 *
 *  Created on: Aug 9, 2013
 *      Author: Henning
 */

#include "ConductanceTree.h"

namespace NetworKit {

Clustering ConductanceTree::bestCutInBfsTree(const Graph& g, node root) {
	count n = g.numberOfNodes();
	std::queue<node> q;
	Clustering visited(n);
	Clustering bestCut(n);
	std::vector<node> parent(n, none);

	edgeweight volume = 0; // volume of set stored in visited
	edgeweight totalVolume = 2 * g.totalEdgeWeight(); // total volume of graph
	edgeweight bndWeight = 0; // weight of cut between visited and remainder of graph
	edgeweight conductance = 0; // conductance of that cut
	edgeweight bestCond = std::numeric_limits<edgeweight>::max(); // best conductance value

	// start with root
	parent[root] = root;
	q.push(root);
	visited.setAll(0);

	auto boundaryWeight([&](edgeweight bndWeight, node v) {
		edgeweight w = bndWeight;

		// subtract degree to vertices within visited
		// add degree to vertices not in visited
		g.forEdgesOf(v, [&](node v, node neighbor) {
			if (visited[neighbor] == 1) {
				w -= g.weight(v, neighbor);
			}
			else {
				w += g.weight(v, neighbor);
			}
		});

		return w;
	});

	auto getConductance([&](edgeweight bndWeight, edgeweight volume) {
		edgeweight denom = std::min(volume, totalVolume - volume);
		return ((double) bndWeight / (double) denom);
	});


	count iterNum = 1;
	while (! q.empty()) {
		node current = q.front();
		q.pop();
		visited[current] = 1;
		volume += g.weightedDegree(current);

		bndWeight = boundaryWeight(bndWeight, current);
		conductance = getConductance(bndWeight, volume);

		DEBUG("conductance: " << conductance);

		if ((conductance < bestCond) && (iterNum < n)) { // do not compute conductance for visited == V
			bestCond = conductance;
			bestCut = visited;
		}

		TRACE("current node in BFS tree: " << current);

		// insert untouched neighbors into queue
		g.forNeighborsOf(current, [&](node neighbor) {
			if (parent[neighbor] == none) {
				q.push(neighbor);
				parent[neighbor] = current;
			}
		});

		++iterNum;
	}

	INFO("Conductance value of best cut found: " << bestCond);

	return bestCut;
}

//node root = 0;
//count n = g.numberOfNodes();
//Clustering clutering = bestCutInBfsTree(g, root);




ConductanceTree::ConductanceTree() {

#if 0
	// create BFS tree with edge weights: conductance value of that cut
	node root = 0;
	BFS bfs;
	std::vector<node> parent = bfs.tree(g, root);

	// check if graph is connected
	g.forNodes([&](node v) {
		if (parent[v] == none) {
			throw std::invalid_argument("Graph is not connected, spanning tree does not exist");
		}
	});

	// if connected: add nodes first
	g.forNodes([&](node v) {
		this->addNode();
	});

	// for each node: add parent edge (except for root)
	g.forNodes([&](node v) {
		this->addEdge(v, parent[v]);
	});
	this->removeEdge(root, root);
#endif
}

ConductanceTree::~ConductanceTree() {

}

} /* namespace NetworKit */
