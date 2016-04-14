/*
 * Matching.cpp
 *
 *  Created on: 03.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "Matching.h"

namespace NetworKit {

Matching::Matching(count z) : data(z, none) {
}

bool Matching::isMatched(const node& u) const {
	return (this->data[u] != none);
}

bool Matching::isProper(const Graph& G) const {
	/**
	 * The content of this data structure represents a matching iff
	 * 	(for all v in V: M[v] = v or M[M[v]] = v) and
	 * 	(for all (u,v) in M): (u,v) in E
	 *
	 */
	bool isProper = true;
	bool sym = true;
	// check if entries are symmetric

	G.forNodes([&](node v) {
		sym = ((data[v] == none) || (data[data[v]] == v));
		if (!sym) {
			DEBUG("node " , v , " is not symmetrically matched");
			isProper = false;
		}
	});

	bool inGraph = true;
	// check if every pair exists as an edge
	G.forNodes([&](node v){
		node w = data[v];
		if ((v != w) && (w != none)) {
			inGraph = G.hasEdge(v, w);
			if (!inGraph) {
				DEBUG("matched pair (" , v , "," , w , ") is not an edge");
				isProper = false;
			}
		}
	});


	return isProper;
}

void Matching::match(const node& u, const node& v) {
	data[u] = v;
	data[v] = u;
}

void Matching::unmatch(const node& u, const node& v) {
	data[u] = none;
	data[v] = none;
}

bool Matching::areMatched(const node& u, const node& v) const {
	return (data[u] == v);
}

count Matching::size(const Graph& G) const {
	count size = 0;
	G.forNodes([&](node v) {
		if (isMatched(v)) {
			++size;
		}
	});
	return size / 2;
}

index Matching::mate(node v) const {
	if (isMatched(v)) {
		return data[v];
	}
	else return none;
}

edgeweight Matching::weight(const Graph& G) const {
	edgeweight weight = 0;

	G.forNodes([&](node v){
		if (isMatched(v) && v < mate(v)) {
			weight += G.weight(v, mate(v));
		}
	});

	return weight;
}

}
 /* namespace NetworKit */
