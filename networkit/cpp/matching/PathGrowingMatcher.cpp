/*
 * PathGrowingMatcher.cpp
 *
 *  Created on: Jun 13, 2013
 *      Author: Henning
 */

#include "PathGrowingMatcher.h"
#include "../auxiliary/PrioQueueForInts.h"

namespace NetworKit {


Matching PathGrowingMatcher::run(Graph& G) {
	// make copy since graph will be transformed
	count n = G.numberOfNodes();

	// init matching to empty
	Matching m1(n);
	Matching m2(n);
	bool takeM1 = true;

	std::vector<count> degrees(n);
	G.forNodes([&](node u) {
		degrees[u] = G.degree(u);
	});
	std::vector<bool> alive(n, true);
	count numEdges = G.numberOfEdges();

//	TRACE("init bpq");
	Aux::PrioQueueForInts bpq(degrees, n-1);

	// main loop
	while (numEdges > 0) {
//		TRACE("Remaining edges: " , numEdges);

		// use vertex with positive degree
		node v = bpq.extractMax();
		assert(v != none);

		// path growing
		while (degrees[v] > 0) {
//			TRACE("Current vertex: " , v);

			// find heaviest incident edge
			node bestNeighbor = 0;
			edgeweight bestWeight = 0;
			G.forEdgesOf(v, [&](node v, node u, edgeweight weight) {
				if (alive[u]) {
					if (weight > bestWeight) {
						bestNeighbor = u;
						bestWeight = weight;
					}
				}
			});

			if (takeM1) {
				// add edge to m1
				m1.match(v, bestNeighbor);
				takeM1 = false;
			}
			else {
				// add edge to m2
				m2.match(v, bestNeighbor);
				takeM1 = true;
			}

			// remove current vertex and its incident edges from graph
//			TRACE("Remove edges of node " , v , ", which has degree " , degrees[v]);

			G.forEdgesOf(v, [&](node v, node u) {
				if (alive[u]) {
					degrees[u]--;
					numEdges--;
					bpq.changePrio(u, degrees[u]);
				}
			});
//			TRACE("Remove node " , v , " of degree " , degrees[v]);

			alive[v] = false;
			bpq.remove(v);

			// start next iteration from best neighbor
			v = bestNeighbor;
		}
	}

	// return the heavier one of the two
	edgeweight weight1 = m1.weight(G);
	edgeweight weight2 = m2.weight(G);
	if (weight1 > weight2)
		return m1;
	else
		return m2;
}

} /* namespace NetworKit */
