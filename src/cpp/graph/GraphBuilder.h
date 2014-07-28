/*
 * GraphBuilder.h
 *
 *  Created on: 15.07.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef GRAPH_BUILDER_H
#define GRAPH_BUILDER_H

#include <vector>

#include "../Globals.h"
#include "Graph.h"

namespace NetworKit {

class GraphBuilder {
protected:
	count n; //!< current number of nodes

	bool weighted; //!< true if the graph will be weighted, false otherwise
	bool directed; //!< true if the graph will be directed, false otherwise

	std::vector< std::vector<node> > halfEdges;
	std::vector< std::vector<edgeweight> > halfWeights;

public:
	GraphBuilder(count n = 0, bool weighted = false, bool directed = false) :
		n(n),
		weighted(weighted),
		directed(directed),
		halfEdges(n),
		halfWeights(weighted ? n : 0)
	{}

	bool isWeighted() const { return weighted; }
	bool isDirected() const { return directed; }
	count numberOfNodes() const { return n; }

	node addNode() {
		halfEdges.push_back(std::vector<node>{});
		if (weighted) {
			halfWeights.push_back(std::vector<edgeweight>{});
		}
		return n++;
	}

	// add edge from u to v
	void addEdge(node u, node v, edgeweight ew = defaultEdgeWeight) {
		halfEdges[u].push_back(v);
		if (weighted) {
			halfWeights[u].push_back(ew);
		}
	}

	/**
	 * Iterate over all undirected pairs of nodes in parallel and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void parallelForNodePairs(L handle) const {
	 	#pragma omp parallel for
		for (node u = 0; u < n; u++) {
			for (node v = u + 1; v < n; v++) {
				handle(u, v);
			}
		}
	}

	Graph toGraph() const {
		Graph G(n, weighted, directed);

		std::vector<count> missingEdgesCounts(n, 0);

		// copy out edges (first half)
		G.forNodes([&](node v) {
			G.outDeg[v] = halfEdges[v].size();
			G.outEdges[v] = halfEdges[v];
			if (weighted) {
				G.outEdgeWeights[v] = halfWeights[v];
			}

			// increase total edge count
			G.m += G.outDeg[v];

			// increase count of incoming edges for neighbors
			for (node u : G.outEdges[v]) {
				missingEdgesCounts[u]++;
			}
		});

		// add second half of the edge
		if (directed) {
			// directed: outEdges is complete, fill inEdges
			// missingEdgesCounts are our inDegrees
			G.inDeg = missingEdgesCounts;

			G.forNodes([&](node v) {
				G.inEdges.reserve(G.inDeg[v]);
			});
			G.forNodes([&](node v) {
				for (node u : halfEdges[v]) {
					// edge v -> u
					G.inEdges[u].push_back(v);
				}
			});

			// same stuff for weights
			if (weighted) {
				G.forNodes([&](node v) {
					G.inEdgeWeights.reserve(G.inDeg[v]);
				});
				G.forNodes([&](node v) {
					for (index i = 0; i < halfWeights[v].size(); i++) {
						node u = halfEdges[v][i];
						edgeweight ew = halfWeights[v][i];
						G.inEdgeWeights[u].push_back(ew);
					}
				});
			}
		} else {
			// undirected: so far each edge is just saved at one node
			// add it to the other node as well
			G.forNodes([&](node v) {
				G.outDeg[v] += missingEdgesCounts[v];
				G.outEdges[v].reserve(G.outDeg[v]);
			});
			G.forNodes([&](node v) {
				for (node u : halfEdges[v]) {
					// edge v -> u
					G.outEdges[u].push_back(v);
				}
			});

			// same stuff for weights
			if (weighted) {
				G.forNodes([&](node v) {
					G.outEdgeWeights[v].reserve(G.outDeg[v]);
				});
				G.forNodes([&](node v) {
					for (index i = 0; i < halfWeights[v].size(); i++) {
						node u = halfEdges[v][i];
						edgeweight ew = halfWeights[v][i];
						G.outEdgeWeights[u].push_back(ew);
					}
				});				
			}
		}

		G.shrinkToFit();
		return G;
	}
};


} /* namespace NetworKit */

#endif /* GRAPH_BUILDER_H */
