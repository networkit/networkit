/*
 * ParallelMatcher.cpp
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "LocalMaxMatcher.h"

namespace NetworKit {

LocalMaxMatcher::LocalMaxMatcher(Graph& graph): Matcher(graph)
{

}

// TODO: update to new edge attribute system

Matching LocalMaxMatcher::run() {
	int64_t z = G.upperNodeIdBound();
	count E = G.numberOfEdges();
	Matching M(z);

	// put edges into array of triples
	struct MyEdge {
		node s;
		node t;
		edgeweight w;
	};

	std::vector<MyEdge> edges(E);
	index e = 0;
	G.forEdges([&](node u, node v, edgeweight w) {
		edges[e].s = u;
		edges[e].t = v;
		edges[e].w = w + Aux::Random::real(1e-3);
		++e;
	});

	// candidates records mating candidates
	std::vector<MyEdge> candidates(z);
	G.parallelForNodes([&](node u) {
		candidates[u].w = (edgeweight) 0;
		candidates[u].t = u; // itself as mating partner => unmatched
	});

	while (E > 0) {
		// for each edge find out if it is locally maximum
		for (auto edge: edges) {
			if (edge.w > candidates[edge.s].w && edge.w > candidates[edge.t].w) {
				candidates[edge.s].s = edge.s;
				candidates[edge.s].t = edge.t;
				candidates[edge.s].w = edge.w;
				candidates[edge.t].s = edge.t;
				candidates[edge.t].t = edge.s;
				candidates[edge.t].w = edge.w;
			}
		}

		// check if candidates agree to match; if so, then match them
		for (auto edge: candidates) {
			node u = edge.s;
			node partner = candidates[u].t;
			if (u < partner && candidates[partner].t == u) {
				// both nodes agree
				M.match(u, partner);
			}
		}


		// create remaining "graph" by selecting remaining edges (as triples)
		// adjust candidates
		std::vector<MyEdge> newEdges;
		for (auto edge: edges) {
			if (! M.isMatched(edge.s) && ! M.isMatched(edge.t)) {
				newEdges.push_back(edge);
				candidates[edge.s].w = (edgeweight) 0;
				candidates[edge.t].w = (edgeweight) 0;
			}
		}
		edges = newEdges;
		E = edges.size();
	}


#if 0
	// local max algorithm
	G.forEdges([&](node u, node v) { // TODO: parallel
		// check neighborhood if both vertices are unmatched
		if (! M.isMatched(u) && ! M.isMatched(v)) {
			edgeweight escore = 0.0;
			bool localMax = true;

			G.forEdgesOf(u, [&](node u, node x) {
				edgeweight otherScore = 0.0;
				if (otherScore > escore) {
					localMax = false;
					//					break;
				}
			});

			G.forEdgesOf(v, [&](node v, node x) {
				edgeweight otherScore = 0.0;

				if (otherScore > escore) {
					localMax = false;
					//					break;
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
#endif

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

		// if nodes mutually prefer each other:
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
