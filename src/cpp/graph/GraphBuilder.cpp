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

GraphBuilder::GraphBuilder(count n, bool weighted, bool directed, bool directSwap) :
	n(n),
	selfloops(0),
	weighted(weighted),
	directed(directed),
	usedirectswap(directSwap),
	halfEdges(n),
	halfEdgeWeights(weighted ? n : 0)
{
	if (directed && directSwap) throw std::runtime_error("Cannot use direct swap in directed graph.");
}

index GraphBuilder::indexHalfEdgeArray(node u, node v) const {
	for (index i = 0; i < halfEdges[u].size(); i++) {
		node x = halfEdges[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

node GraphBuilder::addNode() {
	halfEdges.push_back(std::vector<node>{});
	if (weighted) {
		halfEdgeWeights.push_back(std::vector<edgeweight>{});
	}
	return n++;
}

void GraphBuilder::addEdge(node u, node v, edgeweight ew) {
	halfEdges[u].push_back(v);
	if (weighted) {
		halfEdgeWeights[u].push_back(ew);
	}
	if (u == v) {
	#pragma omp atomic
		selfloops++;
	}
}

void GraphBuilder::setWeight(node u, node v, edgeweight ew) {
	if (!weighted) {
		throw std::runtime_error("Cannot set edge weight in unweighted graph.");
	}

	index vi = indexHalfEdgeArray(u, v);

	if (vi != none) {
		halfEdgeWeights[u][vi] = ew;
	} else {
		addEdge(u, v, ew);
	}
}

void GraphBuilder::increaseWeight(node u, node v, edgeweight ew) {
	if (!weighted) {
		throw std::runtime_error("Cannot increase edge weight in unweighted graph.");
	}

	index vi = indexHalfEdgeArray(u, v);

	if (vi != none) {
		halfEdgeWeights[u][vi] += ew;
	} else {
		addEdge(u, v, ew);
	}
}

Graph GraphBuilder::directSwap() {
	if (directed) throw std::runtime_error("Cannot swap directly in directed Graph.");
	Graph G(n, weighted, directed);
	G.outEdges.swap(halfEdges);
	G.outEdgeWeights.swap(halfEdgeWeights);
	#pragma omp parallel for
	for (node v = 0; v < n; v++) {
		G.outDeg[v] = G.outEdges[v].size();
	}
	correctNumberOfEdges(G, selfloops);

	reset();
	return G;
}

Graph GraphBuilder::toGraphParallel() {
	Graph G(n, weighted, directed);

	// basic idea of the parallelization:
	// 1) each threads collects its own data
	// 2) each node collects all its data from all threads

	int maxThreads = omp_get_max_threads();

	using adjacencylists = std::vector< std::vector<node> >;
	using weightlists = std::vector< std::vector<edgeweight> >;

	std::vector<adjacencylists> inEdgesPerThread(maxThreads, adjacencylists(n));
	std::vector<weightlists> inWeightsPerThread(weighted ? maxThreads : 0, weightlists(n));
	std::vector<count> numberOfSelfLoopsPerThread(maxThreads, 0);

	// step 1
	parallelForNodes([&](node v) {
		int tid = omp_get_thread_num();
		for (index i = 0; i < halfEdges[v].size(); i++) {
			node u = halfEdges[v][i];
			if (directed || u != v) { // self loops don't need to be added twice in undirected graphs
				inEdgesPerThread[tid][u].push_back(v);
				if (weighted) {
					edgeweight ew = halfEdgeWeights[v][i];
					inWeightsPerThread[tid][u].push_back(ew);
				}
			} else {
				numberOfSelfLoopsPerThread[tid]++;
			}
		}
	});

	// we have already half of the edges
	G.outEdges.swap(halfEdges);
	G.outEdgeWeights.swap(halfEdgeWeights);

	// step 2
	parallelForNodes([&](node v) {
		count inDeg = 0;
		count outDeg = G.outEdges[v].size();
		for (int tid = 0; tid < maxThreads; tid++) {
			inDeg += inEdgesPerThread[tid][v].size();
		}

		// allocate enough memory for all edges/weights
		if (directed) {
			G.inEdges[v].reserve(inDeg);
			if (weighted) {
				G.inEdgeWeights[v].reserve(inDeg);
			}
		} else {
			G.outEdges[v].reserve(outDeg + inDeg);
			if (weighted) {
				G.outEdgeWeights[v].reserve(outDeg + inDeg);
			}
		}

		// collect 'second' half of the edges
		if (directed) {
			G.inDeg[v] = inDeg;
			G.outDeg[v] = outDeg;
			for (int tid = 0; tid < maxThreads; tid++) {
				copyAndClear(inEdgesPerThread[tid][v], G.inEdges[v]);
			}
			if (weighted) {
				for (int tid = 0; tid < maxThreads; tid++) {
					copyAndClear(inWeightsPerThread[tid][v], G.inEdgeWeights[v]);
				}	
			}
		} else {
			G.outDeg[v] = inDeg + outDeg;
			for (int tid = 0; tid < maxThreads; tid++) {
				copyAndClear(inEdgesPerThread[tid][v], G.outEdges[v]);
			}
			if (weighted) {
				for (int tid = 0; tid < maxThreads; tid++) {
					copyAndClear(inWeightsPerThread[tid][v], G.outEdgeWeights[v]);
				}	
			}
		}
	});

	// calculate correct m
	count numberOfSelfLoops = 0;
	for (auto c : numberOfSelfLoopsPerThread) {
		numberOfSelfLoops += c;
	}
	correctNumberOfEdges(G, numberOfSelfLoops);

	// bring the builder into an empty, but valid state
	reset();

	return G;
}

Graph GraphBuilder::toGraphSequential() {
	Graph G(n, weighted, directed);

	std::vector<count> missingEdgesCounts(n, 0);
	count numberOfSelfLoops = 0;

	// 'first' half of the edges
	G.outEdges.swap(halfEdges);
	G.outEdgeWeights.swap(halfEdgeWeights);
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
	selfloops = 0;
	halfEdges.clear();
	halfEdgeWeights.clear();
}

void GraphBuilder::correctNumberOfEdges(Graph& G, count numberOfSelfLoops) {
	count edges = 0;
	#pragma omp parallel for reduction(+:edges)
	for (node v = 0; v < G.z; v++) {
		edges += G.degree(v);
	}
	G.m = edges;
	if (!G.isDirected()) {
		// self loops are just counted once
		G.m = numberOfSelfLoops + (G.m - numberOfSelfLoops) / 2;
	}
}

} /* namespace NetworKit */
