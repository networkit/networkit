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
	// TODO: does it pay to parallelize this?
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
	(*this)[u] = v;
	(*this)[v] = u;
}

node& Matching::operator[](const node& u) {
	return this->array[u];
}

const node& Matching::operator [](const node& u) const {
	return this->array[u];
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
