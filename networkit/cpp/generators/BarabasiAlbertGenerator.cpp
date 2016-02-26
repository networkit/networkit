/*
 * BarabasiAlbertGenerator.cpp
 *
 *  Created on: May 28, 2013
 *      Author: forigem
 */

#include "../auxiliary/Random.h"

#include "BarabasiAlbertGenerator.h"

#include <unordered_set>


namespace NetworKit {

BarabasiAlbertGenerator::BarabasiAlbertGenerator() {
}


BarabasiAlbertGenerator::BarabasiAlbertGenerator(count k, count nMax, count n0, bool batagelj) : k(k), nMax(nMax), batagelj(batagelj) {
	if (batagelj) {
		this->n0 = n0;
	} else {
		if (n0 == 0) {
			this->n0 = k;
		}
	}
}

Graph BarabasiAlbertGenerator::generate() {
	if (batagelj) {
		return generateBatagelj();
	}
	Graph G = initializeGraph();
	assert (G.numberOfNodes() >= k);

	for (count i = n0; i < nMax; i++) {
		count degreeSum = G.numberOfEdges() * 2;
		node u = G.addNode();
		std::set<node> targets;
		targets.insert(u);
		int j = 0;
		while (targets.size() - 1 < k) {
			uint64_t random = (uint64_t) Aux::Random::integer(degreeSum);
			j++;
			///if (j > k) throw std::runtime_error("Possible infinite loop detected.");
			bool found = false; // break from node iteration when done
			auto notFound = [&](){ return ! found; };

			G.forNodesWhile(notFound, [&](node v) {


				if (random <= G.degree(v)) {
					found = true; // found a node to connect to
					targets.insert(v);
				}
				random -= G.degree(v);
				//if (j >= G.numberOfNodes() && found==false) throw std::runtime_error("Last node, but still nothing happened.");
			});
		}

		targets.erase(u);

		for (node x : targets) {
			G.addEdge(u, x);
		}

	}

	G.shrinkToFit();
	return G;
}

Graph BarabasiAlbertGenerator::generateBatagelj() {
	count n = nMax;
	Graph G(nMax);
	std::vector<node> M(2 * k * n);
	std::set<std::pair<node, node>> uniqueEdges;
	//std::unordered_set<node> uniqueEdges;

	// initialize n0 connected nodes
	for (index v = 0; v < n0; ++v) {
		M[2 * v ] = v;
		M[2 * v + 1] = v + 1;
	}

	// "draw" the edges
	for (index v = n0; v < n; ++v) {
		for (index i = 0; i < k; ++i) {
			M[2 * (v * k + i)] = v;
			index r = Aux::Random::integer(2 * (v * k + i));
			M[2 * (v * k + i) + 1] = M[r];
		}
	}

	// remove duplicates;
	// for some reason, it is faster with a separate loop when compared to integrating it in the loop aboce
	for (index i = 0; i < (k*n); ++i) {
		//if (!G.hasEdge(M[2*i], M[2*i+1])) G.addEdge(M[2*i], M[2*i+1]);
		uniqueEdges.insert({M[2 * i], M[2 * i + 1]});
		//uniqueEdges.insert(M[2 * i] * nMax + M[2 * i + 1]);
	}
	// add the edges to the graph
	for (const auto& edge : uniqueEdges) {
		G.addEdge(edge.first, edge.second);
		//G.addEdge(edge / nMax, edge % nMax);
	}
	G.shrinkToFit();
	return G;
}

Graph BarabasiAlbertGenerator::initializeGraph() {
	Graph G(0);

	// initialize the graph with n0 connected nodes
	for (count i = 0; i < n0; i++) {
		node v = G.addNode();
		if (i > 0) G.addEdge(v, v - 1);
	}

	return G;
}


} /* namespace NetworKit */
