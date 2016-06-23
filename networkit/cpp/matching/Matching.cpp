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
	return (data.at(u) == v); // TODO: why not also data[v] == u ???
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
	return data.at(v);
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

Partition Matching::toPartition(const Graph& G) const {
	Partition partition(G.upperNodeIdBound());
	std::vector<bool> visited(G.upperNodeIdBound(), false);
	G.forNodes([&](node u){
		if (!visited[u]) {
			if (mate(u) == none) {
				partition.addToSubset(u,u);
			} else {
				partition.addToSubset(u,u);
				partition.addToSubset(u, mate(u));
				visited[u] = true;
				visited[mate(u)] = true;
			}
		}
	});
	return partition;
}

std::vector<node> Matching::getVector() const {
	return this->data; //FIXME is this appropriate? - why not?
}

}
 /* namespace NetworKit */
