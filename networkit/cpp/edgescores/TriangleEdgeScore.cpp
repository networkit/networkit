/*
 * TriangleEdgeScore.cpp
 *
 *  Created on: 29.08.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#include "TriangleEdgeScore.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Timer.h"
#include <omp.h>

namespace NetworKit {

TriangleEdgeScore::TriangleEdgeScore(const Graph& G) : EdgeScore<count>(G) {
}

void TriangleEdgeScore::run() {
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	// direct edge from high to low-degree nodes
	auto isOutEdge = [&](node u, node v) {
		return G.degree(u) > G.degree(v) || (G.degree(u) == G.degree(v) && u < v);
	};

	// Store in-edges explicitly. Idea: all nodes have (relatively) low in-degree
	std::vector<index> inBegin(G.upperNodeIdBound() + 1);
	std::vector<node> inEdges(G.numberOfEdges());

	{
		index pos = 0;
		for (index u = 0; u < G.upperNodeIdBound(); ++u) {
			inBegin[u] = pos;
			if (G.hasNode(u)) {
				G.forEdgesOf(u, [&](node, node v, edgeid eid) {
					if (isOutEdge(v, u)) {
						inEdges[pos++] = v;
					}
				});
			}
		}
		inBegin[G.upperNodeIdBound()] = pos;
	}

	//Edge attribute: triangle count
	std::vector<count> triangleCount(G.upperEdgeIdBound(), 0);
	// Store the id of the out-edges of the current node for fast access later.
	// Due to the parallel iteration this is needed for all thread separately
	// Special values: none = no neighbor. none-1: incoming neighbor.
	std::vector<std::vector<edgeid> > outEdgeId(omp_get_max_threads(), std::vector<count>(G.upperNodeIdBound(), none));

	G.balancedParallelForNodes([&](node u) {
		auto tid = omp_get_thread_num();

		// mark nodes as outgoing neighbor (eid) or incoming neighbor (none-1)
		G.forEdgesOf(u, [&](node _u, node v, edgeid eid) {
			if (isOutEdge(u, v)) {
				outEdgeId[tid][v] = eid;
			} else {
				outEdgeId[tid][v] = none - 1;
			}
		});

		// Find all triangles of the form u-v-w-u where (v, w) is an in-edge.
		// Note that we find each triangle u is part of once.
		G.forEdgesOf(u, [&](node, node v, edgeid eid) {
			// if the edge (u, v) is an out edge (it cannot be none as v is a neighbor)
			bool uvIsOut = (outEdgeId[tid][v] != none - 1);

			// for all in-edges (v, w).
			for (index i = inBegin[v]; i < inBegin[v + 1]; ++i) {
				auto w = inEdges[i];

				edgeid uwid = outEdgeId[tid][w];

				if (uwid != none) {
					// we have found a triangle u-v-w-u

					// record triangle count only for out-edges from u
					// This means that for (v, w) we never record a triangle.
					// The count on (v, w) is updated when w is the central node
					// Record triangle for (u, v)
					if (uvIsOut) {
						++triangleCount[eid];
					}

					// Record triangle for (u, w) if (u, w) is an out-edge, i.e. the out edge id is not none-1
					if (uwid != none - 1) {
						++triangleCount[uwid];
					}
				}
			}
		});

		// Unset all out edge ids
		G.forEdgesOf(u, [&](node, node v, edgeid eid) {
			outEdgeId[tid][v] = none;
		});
	});

	scoreData = std::move(triangleCount);
	hasRun = true;
}

count TriangleEdgeScore::score(node u, node v) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

count TriangleEdgeScore::score(edgeid eid) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
