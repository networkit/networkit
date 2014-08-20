/*
 * GraphBuilder.cpp
 *
 *  Created on: 15.07.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com)
 */

#include <stdexcept>
#include <omp.h>

#include "GraphBuilder.h"

namespace NetworKit {

template <EdgeType edgeType>
GraphBuilder::GraphBuilder(count n, bool weighted, bool directed) :
	n(n),
	weighted(weighted),
	directed(directed),
	outEdges(n),
	inEdges(edgeType == GraphBuilderEdgeType::FullEdges ? n : 0),
	outEdgeWeights(weighted ? n : 0),
	inEdgeWeights(weighted && edgeType == GraphBuilderEdgeType::FullEdges ? n : 0) {
}

template <EdgeType edgeType>
index GraphBuilder::indexInOutEdgeArray(node u, node v) const {
	for (index i = 0; i < outEdges[u].size(); i++) {
		node x = outEdges[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

template <EdgeType edgeType>
node GraphBuilder::addNode() {
	outEdges.push_back(std::vector<node>{});
	if (weighted) {
		outEdgeWeights.push_back(std::vector<edgeweight>{});
	}
	return n++;
}

template <EdgeType edgeType>
void GraphBuilder::addEdge(node u, node v, edgeweight ew) {
	outEdges[u].push_back(v);
	if (weighted) {
		outEdgeWeights[u].push_back(ew);
	}
}

template<>
void GraphBuilder<EdgeType::HalfEdges>::addInEdge(node u, node v, edgeweight) = delete;

template<>
void GraphBuilder<EdgeType::FullEdges>::addInEdge(node u, node v, edgeweight) {
	if (!directed()) {
		throw std::runtime_error("Cannot set in-edges in undirected graphs, use addEdge.");
	}
	inEdges[v].push_back(u);
	if (weighted) {
		inEdgeWeights[v].push_back(u);
	}
}

void GraphBuilder::setWeight(node u, node v, edgeweight ew) {
	if (!weighted) {
		throw std::runtime_error("Cannot set edge weight in unweighted graph.");
	}

	index vi = indexInOutEdgeArray(u, v);

	if (vi != none) {
		outEdgeWeights[u][vi] = ew;
	} else {
		addEdge(u, v, ew);
	}
}

void GraphBuilder::increaseWeight(node u, node v, edgeweight ew) {
	if (!weighted) {
		throw std::runtime_error("Cannot increase edge weight in unweighted graph.");
	}

	index vi = indexInOutEdgeArray(u, v);

	if (vi != none) {
		outEdgeWeights[u][vi] += ew;
	} else {
		addEdge(u, v, ew);
	}
}

template <GraphBuilderEdgeType edgeType>
Graph GraphBuilder::toGraph_directSwap() {
	Graph G(n, weighted, directed);
	G.outEdges.swap(outEdges);
	G.outEdges.swap(inEdges);
	G.outEdgeWeights.swap(outEdgeWeights);
	G.outEdgeWeights.swap(inEdgeWeights);
	
	count selfLoopCount = 0;

	#pragma omp parallel for reduction(+:selfLoopCount)
	for (node v = 0; v < n; v++) {
		G.outDeg[v] = G.outEdges[v].size();
		if (directed) {
			G.inDeg[v] = G.inEdges[v].size();
		}
		selfLoopCount += std::count(G.outEdges[v].begin(), G.outEdges[v].end(), v);
	}
	correctNumberOfEdges(G, selfLoopCount);

	reset();
	return G;
}

template <GraphBuilderEdgeType edgeType>
Graph GraphBuilder::toGraph_sequential() {
	Graph G(n, weighted, directed);

	std::vector<count> missingEdgesCounts(n, 0);
	count numberOfSelfLoops = 0;

	// 'first' half of the edges
	G.outEdges.swap(outEdges);
	G.outEdgeWeights.swap(outEdgeWeights);
	parallelForNodes([&](node v) {
		G.outDeg[v] = G.outEdges[v].size();
	});

	// count missing edges for each node
	G.forNodes([&](node v) {
		// increase count of incoming edges for all neighbors
		for (node u : G.outEdges[v]) {
			if (directed || u != v) {
				missingEdgesCounts[u]++;
			} else {
				// self loops don't need to be added again
				// but we need to count them to correct the number of edges later
				numberOfSelfLoops++;
			}
		}
	});

	// 'second' half the edges
	if (directed) {
		// directed: outEdges is complete, missing half edges are the inEdges
		// missingEdgesCounts are our inDegrees
		G.inDeg = missingEdgesCounts;

		// reserve the exact amount of space needed
		G.forNodes([&](node v) {
			G.inEdges[v].reserve(G.inDeg[v]);
			if (weighted) {
				G.inEdgeWeights[v].reserve(G.inDeg[v]);
			}
		});

		// copy values
		G.forNodes([&](node v) {
			for (index i = 0; i < G.outDeg[v]; i++) {
				node u = G.outEdges[v][i];
				G.inEdges[u].push_back(v);
				if (weighted) {
					edgeweight ew = G.outEdgeWeights[v][i];
					G.inEdgeWeights[u].push_back(ew);
				}
			}
		});
	} else {
		// undirected: so far each edge is just saved at one node
		// add it to the other node as well

		// reserve the exact amount of space needed
		G.forNodes([&](node v) {
			G.outEdges[v].reserve(G.outDeg[v] + missingEdgesCounts[v]);
			if (weighted) {
				G.outEdgeWeights[v].reserve(G.outDeg[v] + missingEdgesCounts[v]);
			}
		});

		// cope value
		G.forNodes([&](node v) {
			// the first G.outDeg[v] edges in G.outEdges[v] are the first half edges
			// we are adding after G.outDeg[v]
			for (index i = 0; i < G.outDeg[v]; i++) {
				node u = G.outEdges[v][i];
				if (u != v) {
					G.outEdges[u].push_back(v);
					if (weighted) {
						edgeweight ew = G.outEdgeWeights[v][i];
						G.outEdgeWeights[u].push_back(ew);
					}
				} else {
					// ignore self loops here
				}
			}
		});

		// correct degree
		G.forNodes([&](node v) {
			G.outDeg[v] += missingEdgesCounts[v];
		});
	}

	// calculate correct m	
	correctNumberOfEdges(G, numberOfSelfLoops);

	// bring the builder into an empty, but valid state
	reset();

	return G;
}

void GraphBuilder::reset() {
	n = 0;
	outEdges.clear();
	inEdges.clear();
	outEdgeWeights.clear();
	inEdgeWeights.clear();
}

void GraphBuilder::correctNumberOfEdges(Graph& G, count numberOfSelfLoops) {
	count edgeCount = 0;
	#pragma omp parallel for reduction(+:edgeCount)
	for (node v = 0; v < G.z; v++) {
		edgeCount += G.degree(v);
	}
	G.m = edgeCount;
	if (!G.isDirected()) {
		// self loops are just counted once
		G.m = numberOfSelfLoops + (G.m - numberOfSelfLoops) / 2;
	}
}

} /* namespace NetworKit */
