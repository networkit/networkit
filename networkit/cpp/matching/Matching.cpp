/*
 * Matching.cpp
 *
 *  Created on: 03.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "Matching.h"

namespace NetworKit {

Matching::Matching(uint64_t n) : data(n, none), n(n) {
}

bool Matching::isMatched(const node& u) const {
	return (this->data[u] != none);
}

bool Matching::isProper(Graph& G) const {
	/**
	 * The content of this data structure represents a matching iff
	 * 	(for all v in V: M[v] = v or M[M[v]] = v) and
	 * 	(for all (u,v) in M): (u,v) in E
	 *
	 */
	bool sym = true;
	// check if entries are symmetric
	for (node v = 0; v < G.numberOfNodes(); ++v) {
		sym = ((data[v] == none) || (data[data[v]] == v));
		if (!sym) {
			DEBUG("node " , v , " is not symmetrically matched");
			return false;
		}
	}

	bool inGraph = true;
	// check if every pair exists as an edge
	for (node v = 0; v < G.numberOfNodes(); ++v) {
		node w = data[v];
		if ((v != w) && (w != none)) {
			inGraph = G.hasEdge(v, w);
			if (!inGraph) {
				DEBUG("matched pair (" , v , "," , w , ") is not an edge");
				return false;
			}
		}
	}

	return (sym && inGraph);
}

void Matching::match(const node& u, const node& v) {
	data[u] = v;
	data[v] = u;
}

void Matching::unmatch(const node& u, const node& v) {
	data[u] = u;
	data[v] = v;
}

bool Matching::areMatched(const node& u, const node& v) const {
	return (data[u] == v);
}

count Matching::size() const {
	count size = 0;
	for (index i = 0; i < n; ++i) { // TODO: parallel
		if (isMatched(i)) {
			++size;
		}
	}
	return size / 2;
}

index Matching::mate(node v) const {
	if (isMatched(v)) {
		return data[v];
	}
	else return none;
}

edgeweight Matching::weight(const Graph& g) const {
	edgeweight weight = 0;

	for (index i = 0; i < n; ++i) {
		if (isMatched(i) && i < mate(i)) {
			weight += g.weight(i, mate(i));
		}
	}

	return weight;
}

}
 /* namespace NetworKit */
