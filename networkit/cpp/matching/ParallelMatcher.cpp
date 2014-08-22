/*
 * ParallelMatcher.cpp
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ParallelMatcher.h"

namespace NetworKit {

LocalMaxMatcher::LocalMaxMatcher(uint64_t attrId_) :
		attrId(attrId_) {

}

Matching LocalMaxMatcher::run(Graph& graph) {
	Graph G = graph;

	int64_t n = G.numberOfNodes();
	Matching M(n);

	// local max algorithm
	G.forEdges([&](node u, node v) { // TODO: parallel
		// check neighborhood if both vertices are unmatched
		if (! M.isMatched(u) && ! M.isMatched(v)) {
			edgeweight escore = 0.0;
			// FIXME: update to new edge attribute system
			// edgeweight escore = (attrId == none) ?
			// 		(G.weight(u, v)):
			// 		(G.attribute_double(u, v, attrId));
			bool localMax = true;

			G.forEdgesOf(u, [&](node u, node x) {
				edgeweight otherScore = 0.0;
				// FIXME: update to new edge attribute system
				// edgeweight otherScore = (attrId == none) ?
				// 		(G.weight(u, x)):
				// 		(G.attribute_double(u, x, attrId));

				if (otherScore > escore) {
					localMax = false;
// TODO					break;
				}
			});

			G.forEdgesOf(v, [&](node v, node x) {
				edgeweight otherScore = 0.0;
				// FIXME: update to new edge attribute system
				// edgeweight otherScore = (attrId == none) ?
				// 		(G.weight(v, x)):
				// 		(G.attribute_double(v, x, attrId));

				if (otherScore > escore) {
					localMax = false;
// TODO					break;
				}
			});

			if (localMax) {
				M.match(u, v);

				// remove incident edges
				G.forEdgesOf(u, [&](node u, node x) {
					G.removeEdge(u, x);
				});

				G.forEdgesOf(v, [&](node v, node x) {
					G.removeEdge(v, x);
				});
			}
		}
	});

#if 0
	// TODO: exclude isolated nodes?

	// FIXME: will crash if deleted nodes present
	int64_t n = G.numberOfNodes();
	std::vector<node> candidate(n, 0);//!< candidate[v] is the preferred matching partner of v
	std::vector<std::set<node> > S(n, std::set<node>());//!< S[v] is a set with the potential
													//!< candidates of node v
	std::vector<node> D;							//!< targets of dominating edges
	Matching M(n);

	G.forNodes([&](node v) {
				// S[v] <- N(v)
				G.forNeighborsOf(v, [&](node w) {
							S[v].insert(w);
						});
				// set candidate of v to neighbor x with strongest connection
				// INFO: argmax is equivalent to Python max(collection, key=function)
				auto cv = argmax(S[v], [&](node w) {
							return G.weight(v, w);
						});
				candidate[v] = cv;

				// if nodes mutually prefer eachother:
				if (candidate[candidate[v]] == v) {
					D.push_back(v);
					D.push_back(candidate[v]);
					M.match(v, candidate[v]);
				}
			});

	while (! D.empty()) {
		node v = D.back();
		D.pop_back();

		for (node x : S[v]) {
			if ((x != candidate[v]) && (! M.areMatched(x, candidate[x]))) {
				S[x].erase(v);
				if (! S[x].empty()) {
					// find new candidate
					auto cx = argmax(S[x], [&](node y) {
								return G.weight(x, y);
							});
					candidate[x] = cx;
				} else {
					// TODO: what if no potential candidates are left?
				}
				if (candidate[candidate[x]] == x) {
					D.push_back(x);
					D.push_back(candidate[x]);
					M.match(x, candidate[x]);
				}

			}
		}
	} // end while
#endif

	return M;
}

} /* namespace NetworKit */
