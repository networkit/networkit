/*
 * Matching.cpp
 *
 *  Created on: 03.12.2012
 *      Author: cls
 */

#include "Matching.h"

namespace EnsembleClustering {

Matching::Matching(int n) {
	this->array = new node[n];
	// initialize each node's matching partner to itself
	for (int i = 0; i < n; ++i) {
		this->array[i] = i;
	}
}

Matching::~Matching() {
	delete[] this->array;
}


bool Matching::isMatched(const node& u) const {
	return (this->array[u] != u);
}

bool Matching::isProper(Graph& G) const {
}



void Matching::match(const node& u, const node& v) {
	this->array[u] = v;
	this->array[v] = u;
}

node& Matching::operator[](const node& u) {
	return this->array[u];
}

const node& Matching::operator [](const node& u) const {
	return this->array[u];
}

} /* namespace EnsembleClustering */
