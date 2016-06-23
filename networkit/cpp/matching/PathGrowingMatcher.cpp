/*
 * PathGrowingMatcher.cpp
 *
 *  Created on: Jun 13, 2013
 *      Author: Henning
 */

#include "PathGrowingMatcher.h"
#include "../auxiliary/PrioQueueForInts.h"

namespace NetworKit {

PathGrowingMatcher::PathGrowingMatcher(const Graph& G): Matcher(G) {
	if (G.numberOfSelfLoops() > 0) {
		throw std::invalid_argument("G has self-loops and cannot be processed");
	}
}


PathGrowingMatcher::PathGrowingMatcher(const Graph& G, const std::vector<double>& edgeScores): Matcher(G, edgeScores) {
	if (G.numberOfSelfLoops() > 0) {
		throw std::invalid_argument("G has self-loops and cannot be processed");
	}
}

void PathGrowingMatcher::run() {
	count z = G.upperNodeIdBound();

	// init matching to empty
	Matching m1(z);
	Matching m2(z);
	bool takeM1 = true;

	// degrees tracks degree of vertices,
	// avoids to make a copy of the graph and
	// delete vertices and edges explicitly.
	// Init to none important because deleted nodes
	// must have value none for priority queue.
	std::vector<count> degrees(z, none);
	G.parallelForNodes([&](node u) {
		degrees.at(u) = G.degree(u);
	});

	// alive tracks if vertices are alive or not in the algorithm
	std::vector<bool> alive(z, true);
	count numEdges = G.numberOfEdges();

	// PQ to retrieve vertices with degree > 0 quickly
	Aux::PrioQueueForInts bpq(degrees, z-1);

	// main loop
	while (numEdges > 0) {
		// use vertex with positive degree
		node v = bpq.extractMax();
		assert(v != none);

		// path growing
		while (degrees[v] > 0) {
			// find heaviest incident edge
			node bestNeighbor = 0;
			edgeweight bestWeight = 0;

			if (edgeScoresAsWeights) {
				G.forEdgesOf(v, [&](node v, node u, edgeid eid) {
					if (alive.at(u)) {
						if (edgeScores.at(eid) > bestWeight) {
							bestNeighbor = u;
							bestWeight = edgeScores.at(eid);
						}
					}
				});
			} else {
				G.forEdgesOf(v, [&](node v, node u, edgeweight weight) {
					if (alive.at(u)) {
						if (weight > bestWeight) {
							bestNeighbor = u;
							bestWeight = weight;
						}
					}
				});
			}


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
			G.forEdgesOf(v, [&](node v, node u) {
				if (alive.at(u)) {
					degrees.at(u)--;
					numEdges--;
					bpq.changePrio(u, degrees.at(u));
				}
			});
			alive[v] = false;
			bpq.remove(v);

			// start next iteration from best neighbor
			v = bestNeighbor;
		}
	}

	// return the heavier one of the two
	edgeweight weight1;
	if (edgeScoresAsWeights) {
		G.forEdges([&](node u, node v, edgeid eid){
			weight1 += edgeScores.at(eid);
		});
	} else {
		weight1 = m1.weight(G);
	}
	edgeweight weight2;
	if (edgeScoresAsWeights) {
		G.forEdges([&](node u, node v, edgeid eid){
			weight2 += edgeScores.at(eid);
		});
	} else {
		weight2 = m2.weight(G);
	}
	INFO("weight of first matching: ", weight1);
	INFO("weight of second matching: ", weight2);

	if (weight1 > weight2)
		M = m1;
	else
		M = m2;

	hasRun = true;
}

} /* namespace NetworKit */
