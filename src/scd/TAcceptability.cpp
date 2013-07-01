/*
 * TAcceptability.cpp
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#include "TAcceptability.h"

namespace NetworKit {

TAcceptability::TAcceptability(const Graph& G, std::unordered_set<node>& community, std::unordered_set<node>& shell) : G(&G), community(&community), shell(&shell) {
	// TODO Auto-generated constructor stub

}

TAcceptability::~TAcceptability() {
	// TODO Auto-generated destructor stub
}

double TAcceptability::getValue(node v) {
	throw std::runtime_error("calling abstract method - needs to be redefined in subclass");
}

TDummyAcceptability::TDummyAcceptability(const Graph& G,
		std::unordered_set<node>& community, std::unordered_set<node>& shell) : TAcceptability(G, community, shell) {
}

TDummyAcceptability::~TDummyAcceptability() {
}

double TDummyAcceptability::getValue(node v) {
	return 0.5;
}

TNodeClusterSimilarity::TNodeClusterSimilarity(
	const Graph& G, std::unordered_set<node>& community, std::unordered_set<node>& shell): TAcceptability(G, community, shell) {
}

TNodeClusterSimilarity::~TNodeClusterSimilarity() {
}

double TNodeClusterSimilarity::getValue(node v) {

	double intersection = 0;
	this->G->forNeighborsOf(v, [&](node u) {

		if (this->community->find(u) != this->community->end()||this->shell->find(u) != this->shell->end()) {
			intersection++;
		}
	});

	if (this->shell->find(v) != this->shell->end()) {
		intersection++;
	}

	if (G->hasEdge(v, v)) {
		return (intersection - 1) / (G->degree(v) + community->size() + shell->size() - intersection  + 1);
	}

	return intersection / (G->degree(v) + 1 + community->size() + shell->size() - intersection);
}
} /* namespace NetworKit */


