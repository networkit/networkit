/*
 * Matching.cpp
 *
 *  Created on: 03.12.2012
 *      Author: cls
 */

#include "Matching.h"

namespace EnsembleClustering {

Matching::Matching(int64_t n) : NodeMap<node>(n, 0) {
	// initialize each node's matching partner to itself
	for (int64_t i = 1; i < n+1; ++i) {
		this->array[i] = i;
	}
}

Matching::~Matching() {
}


bool Matching::isMatched(const node& u) const {
	return (this->array[u] != u);
}

bool Matching::isProper(Graph& G) const {
	/**
	 * The content of this data structure represents a matching iff
	 * 	(for all v in V: M[v] = 0 or M[M[v]] = v) and
	 * 	(for all (u,v) in M): (u,v) in E
	 *
	 */
	bool sym;
	// check if entries are symmetric
	for (node v = 1; v <= G.numberOfNodes(); ++v) {
		sym = (((*this)[v] == 0) || ((*this)[(*this)[v]] == v));
		if (!sym) {
			DEBUG("node " << v << " is not symmetrically matched");
			return false;
		}
	}
	bool inGraph;
	// check if every pair exists as an edge
	for (node v = 1; v <= G.numberOfNodes(); ++v) {
		node w = (*this)[v];
		inGraph = G.hasEdge(v, w);
		if (!inGraph) {
			DEBUG("matched pair (" << v << "," << w << ") is not an edge");
			return false;
		}
	}

}


void Matching::match(const node& u, const node& v) {
	(*this)[u] = v;
	(*this)[v] = u;
}


Matching& Matching::operator =(const Matching& from) {
	// clone and dispose properly on assignment
	if (&from != this) {
		this->dispose();
		this->clone(from);
	}
	return *this;
}

void Matching::clone(const Matching& from) {
	// TODO: deep copy array
}

void Matching::dispose() {
	delete[] this->array;
}

} /* namespace EnsembleClustering */
