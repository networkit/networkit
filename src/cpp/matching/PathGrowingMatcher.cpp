/*
 * PathGrowingMatcher.cpp
 *
 *  Created on: Jun 13, 2013
 *      Author: Henning
 */

#include "PathGrowingMatcher.h"

namespace NetworKit {


Matching PathGrowingMatcher::run(Graph& G) {
	// make copy since graph will be transformed
	// copy graph because we make changes due to merges
	Graph graph(G.numberOfNodes(), true); // make weighted copy
	G.forEdges([&](node u, node v, edgeweight w){
		graph.addEdge(u, v, w);
	});
	const count n = graph.numberOfNodes();

	// init matching to empty
	Matching m1(n);
	Matching m2(n);
	bool smaller = true;

	// main loop
	while (graph.numberOfEdges() > 0) {
		TRACE("Remaining edges: " , graph.numberOfEdges());

		// use arbitrary vertex with positive degree
		node v = 0;
		while (graph.degree(v) == 0 && v < n) {
			++v;
		}

		// path growing
		while (graph.degree(v) > 0) {
			TRACE("Current vertex: " , v);

			// find heaviest incident edge
			node bestNeighbor = 0;
			edgeweight bestWeight = 0;
			graph.forEdgesOf(v, [&](node v, node u, edgeweight weight) {
				if (weight > bestWeight) {
					bestNeighbor = u;
					bestWeight = weight;
				}
			});

			if (smaller) {
				// add edge to m1
				m1.match(v, bestNeighbor);
				smaller = false;
			}
			else {
				// add edge to m2
				m2.match(v, bestNeighbor);
				smaller = true;
			}

			// remove current vertex and its incident edges from graph
			TRACE("Remove edges of node " , v , ", which has degree " , graph.degree(v));
			graph.forEdgesOf(v, [&](node v, node u) {
				graph.removeEdge(v, u);
			});
			TRACE("Remove node " , v , " of degree " , graph.degree(v));
			graph.removeNode(v);

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
