/*
 * BarabasiAlbertGenerator.cpp
 *
 *  Created on: May 28, 2013
 *      Author: forigem
 */

#include "../auxiliary/Random.h"

#include "BarabasiAlbertGenerator.h"

#include <set>


namespace NetworKit {

BarabasiAlbertGenerator::BarabasiAlbertGenerator() {
}


BarabasiAlbertGenerator::BarabasiAlbertGenerator(count k, count nMax, count n0, bool batagelj) : initGraph(0), k(k), nMax(nMax), batagelj(batagelj) {
	if (k > nMax)
		throw std::runtime_error("k (number of attachments per node) may not be larger than the number of nodes in the target graph (nMax)");
	if (n0 > nMax)
		throw std::runtime_error("n0 (number of initially connected nodes) may not be larger than the number of nodes in the target graph (nMax)");
	if (batagelj) {
		this->n0 = n0;
	} else {
		if (n0 < k) {
			if (n0 > 0) {
				WARN("given n0 is smaller than k, setting n0 = k");
			}
			this->n0 = k;
		} else {
			this->n0 = n0;
		}
	}
}

BarabasiAlbertGenerator::BarabasiAlbertGenerator(count k, count nMax, const Graph& initGraph, bool batagelj) : initGraph(initGraph), k(k), nMax(nMax), n0(0), batagelj(batagelj) {
	if (initGraph.numberOfNodes() != initGraph.upperNodeIdBound())
		throw std::runtime_error("initGraph is expected to have consecutive node ids");
	if (k > nMax)
		throw std::runtime_error("k (number of attachments per node) may not be larger than the number of nodes in the target graph (nMax)");
	if (initGraph.numberOfNodes() > nMax)
		throw std::runtime_error("initialization graph cannot have more nodes than the target graph (nMax)");
	if (!batagelj && initGraph.numberOfNodes() < k) {
		throw std::runtime_error("initialization graph for the original method needs at least k nodes");
	}
}

Graph BarabasiAlbertGenerator::generate() {
	if (batagelj) {
		return generateBatagelj();
	} else {
		return generateOriginal();
	}
}

Graph BarabasiAlbertGenerator::generateOriginal() {
	Graph G(nMax);
	if (n0 != 0) {
		// initialize the graph with n0 connected nodes
		for (count i = 1; i < n0; i++) {
			G.addEdge(i-1, i);
		}
	} else {
		// initialize the graph with the edges from initGraph
		// and set n0 accordingly
		initGraph.forEdges([&G](node u, node v) {
			G.addEdge(u, v);
		});
		n0 = initGraph.upperNodeIdBound();
	}
	assert (G.numberOfNodes() >= k);

	for (count i = n0; i < nMax; i++) {
		count degreeSum = G.numberOfEdges() * 2;
		node u = i;
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

	if (initGraph.numberOfNodes() == 0) {
		// initialize n0 connected nodes
		for (index v = 0; v < n0; ++v) {
			M[2 * v ] = v;
			M[2 * v + 1] = v + 1;
		}
	} else {
		index i = 0;
		initGraph.forEdges( [&M,&i] (node u, node v) {
			M[2 * i] = u;
			M[2 * i +1] = v;
			++i;
		});
		n0 = i;
	}

	// "draw" the edges
	for (index v = n0; v < n; ++v) {
		for (index i = 0; i < k; ++i) {
			M[2 * (v * k + i)] = v;
			index r = Aux::Random::integer(2 * (v * k + i));
			M[2 * (v * k + i) + 1] = M[r];
		}
	}

	// remove duplicates and avoid selfloops
	// for some reason, it seems to be faster with a separate loop when compared to integrating it in the loop above
	for (index i = 0; i < (k*n); ++i) {
		if (M[2 * i] != M[2 * i + 1])
			uniqueEdges.insert( std::minmax({M[2 * i], M[2 * i + 1]}) );
	}
	// add the edges to the graph
	for (const auto& edge : uniqueEdges) {
		G.addEdge(edge.first, edge.second);
	}
	G.shrinkToFit();
	return G;
}

} /* namespace NetworKit */
