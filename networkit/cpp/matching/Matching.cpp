/*
 * Matching.cpp
 *
 *  Created on: 03.12.2012
 */

#include "Matching.h"

namespace NetworKit {


Matching::Matching(count z) : data(z, none) {
}

bool Matching::isMatched(node u) const {
	return (this->data.at(u) != none);
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
		sym = ((data.at(v) == none) || (data[data.at(v)] == v));
		if (!sym) {
			DEBUG("node " , v , " is not symmetrically matched");
			isProper = false;
		}
	});

	bool inGraph = true;
	// check if every pair exists as an edge
	G.forNodes([&](node v){
		node w = data.at(v);
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

void Matching::match(node u, node v) {
	data.at(u) = v;
	data.at(v) = u;
}

void Matching::unmatch(node u, node v) {
	data.at(u) = none;
	data.at(v) = none;
}

bool Matching::areMatched(node u, node v) const {
	return (data.at(u) == v);
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
		return data.at(v);
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
