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

GraphBuilder::GraphBuilder(count n, bool weighted, bool directed) :
	n(n),
	weighted(weighted),
	directed(directed),
	halfEdges(n),
	halfEdgeWeights(weighted ? n : 0)
{
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
}

void GraphBuilder::increaseWeight(node u, node v, edgeweight ew) {
	if (!weighted) {
		throw std::runtime_error("Cannot increase edge weight in unweighted graph.");
	}

	index vi = indexHalfEdgeArray(u, v);

	if (directed) {
		if (vi != none) {
			halfEdgeWeights[u][vi] += ew;
		} else {
			addEdge(u, v, ew);
		}
	} else {
		if (vi != none) {
			halfEdgeWeights[u][vi] += ew;
		} else {
			index ui = indexHalfEdgeArray(v, u);
			if (ui != none) {
				halfEdgeWeights[v][ui] += ew;
			} else {
				addEdge(u, v, ew);
			}
		}
	}
}

Graph GraphBuilder::toGraphParallel() {
	Graph G(n, weighted, directed);

	int maxThreads = omp_get_max_threads();

	using adjacencylists = std::vector< std::vector<node> >;
	using weightlists = std::vector< std::vector<edgeweight> >;

	std::vector<adjacencylists> inEdgesPerThread(maxThreads, adjacencylists(n));
	std::vector<weightlists> inWeightsPerThread(weighted ? maxThreads : 0, weightlists(n));

	parallelForNodes([&](node v) {
		int tid = omp_get_thread_num();
		for (node u : halfEdges[v]) {
			inEdgesPerThread[tid][u].push_back(v);
		}
		if (weighted) {
			for (index i = 0; i < halfEdges[v].size(); i++) {
				node u = halfEdges[v][i];
				edgeweight ew = halfEdgeWeights[v][i];
				inWeightsPerThread[tid][u].push_back(ew);
			}
		}
	});

	parallelForNodes([&](node v) {
		count inDeg = 0;
		count outDeg = halfEdges[v].size();
		for (int tid = 0; tid < maxThreads; tid++) {
			inDeg += inEdgesPerThread[tid][v].size();
		}

		if (directed) {
			G.inEdges[v].reserve(inDeg);
			G.outEdges[v].reserve(outDeg);
			if (weighted) {
				G.inEdgeWeights[v].reserve(inDeg);
				G.outEdgeWeights[v].reserve(outDeg);
			}
		} else {
			G.outEdges[v].reserve(outDeg + inDeg);
			if (weighted) {
				G.outEdgeWeights[v].reserve(outDeg + inDeg);
			}
		}

		std::copy(halfEdges[v].begin(), halfEdges[v].end(), G.outEdges[v].end());
		halfEdges[v].clear();
		if (weighted) {
			std::copy(halfEdgeWeights[v].begin(), halfEdgeWeights[v].end(), G.outEdgeWeights[v].end());
			halfEdgeWeights[v].clear();
		}

		if (directed) {
			G.inDeg[v] = inDeg;
			G.outDeg[v] = outDeg;
			for (int tid = 0; tid < maxThreads; tid++) {
				std::copy(inEdgesPerThread[tid][v].begin(), inEdgesPerThread[tid][v].end(), G.inEdges[v].end());
				inEdgesPerThread[tid][v].clear();
			}
			if (weighted) {
				for (int tid = 0; tid < maxThreads; tid++) {
					std::copy(inWeightsPerThread[tid][v].begin(), inWeightsPerThread[tid][v].end(), G.inEdgeWeights[v].end());
					inWeightsPerThread[tid][v].clear();
				}	
			}
		} else {
			G.outDeg[v] = inDeg + outDeg;
			for (int tid = 0; tid < maxThreads; tid++) {
				std::copy(inEdgesPerThread[tid][v].begin(), inEdgesPerThread[tid][v].end(), G.outEdges[v].end());
				inEdgesPerThread[tid][v].clear();
			}
			if (weighted) {
				for (int tid = 0; tid < maxThreads; tid++) {
					std::copy(inWeightsPerThread[tid][v].begin(), inWeightsPerThread[tid][v].end(), G.outEdgeWeights[v].end());
					inWeightsPerThread[tid][v].clear();
				}	
			}
		}
	});

	// calculate correct m
	forNodes([&](node v) {
		G.m += G.degree(v);
	});
	if (!directed) {
		G.m /= 2;
	}

	// bring the builder into an empty, but valid state
	n = 0;
	halfEdges.clear();
	halfEdgeWeights.clear();

	return G;
}

Graph GraphBuilder::toGraphSequential() {
	Graph G(n, weighted, directed);

	std::vector<count> missingEdgesCounts(n, 0);

	// first half edge
	// copy halfEdges to G.outEdges and set G.outDeg
	G.forNodes([&](node v) {
		G.outDeg[v] = halfEdges[v].size();
		G.outEdges[v] = halfEdges[v];
		halfEdges[v].clear();
	});

	// same for weights
	if (weighted) {
		G.forNodes([&](node v) {
			G.outEdgeWeights[v] = halfEdgeWeights[v];
			halfEdgeWeights[v].clear();
		});		
	}

	// count missing edges for each node
	G.forNodes([&](node v) {
		// increase count of incoming edges for all neighbors
		for (node u : G.outEdges[v]) {
			missingEdgesCounts[u]++;
		}
	});

	// second half edge
	if (directed) {
		// directed: outEdges is complete, missing half edges are the inEdges
		// missingEdgesCounts are our inDegrees
		G.inDeg = missingEdgesCounts;

		// reserve the exact amount of space needed first
		G.forNodes([&](node v) {
			G.inEdges[v].reserve(G.inDeg[v]);
			if (weighted) {
				G.inEdgeWeights[v].reserve(G.inDeg[v]);
			}
		});
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

		// reserve the exact amount of space needed first
		G.forNodes([&](node v) {
			G.outEdges[v].reserve(G.outDeg[v] + missingEdgesCounts[v]);
			if (weighted) {
				G.outEdgeWeights[v].reserve(G.outDeg[v] + missingEdgesCounts[v]);
			}
		});
		G.forNodes([&](node v) {
			// the first G.outDeg[v] edges in G.outEdges[v] are the first half edges
			// we are adding after G.outDeg[v]
			for (index i = 0; i < G.outDeg[v]; i++) {
				node u = G.outEdges[v][i];
				G.outEdges[u].push_back(v);
				if (weighted) {
					edgeweight ew = G.outEdgeWeights[v][i];
					G.outEdgeWeights[u].push_back(ew);
				}
			}
		});

		// correct degree
		G.forNodes([&](node v) {
			G.outDeg[v] += missingEdgesCounts[v];
		});
	}

	// calculate correct m	
	G.forNodes([&](node v) {
		G.m += G.degree(v);
	});
	if (!directed) {
		G.m /= 2;
	}

	// bring the builder into an empty, but valid state
	n = 0;
	halfEdges.clear();
	halfEdgeWeights.clear();

	return G;
}

} /* namespace NetworKit */
