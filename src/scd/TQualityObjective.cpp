/*
 * TQualityObjective.cpp
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#include "TQualityObjective.h"

namespace NetworKit {

} /* namespace NetworKit */

NetworKit::TQualityObjective::TQualityObjective(const Graph& G,
		std::unordered_set<node>& community) : G(&G), community(&community), volume(0), nBoundaryEdges(0) {
	this->degSum = this->G->parallelSumForNodes([&](node u){
		return this->G->degree(u);
	});
}

NetworKit::TQualityObjective::~TQualityObjective() {
}

double NetworKit::TQualityObjective::getValue(node v) {
	throw std::runtime_error("abstract method called - must be redefined in subclass");
}

NetworKit::TLocalModularityM::TLocalModularityM(const Graph& G,
		std::unordered_set<node>& community) : TQualityObjective(G, community) {
}

NetworKit::TLocalModularityM::~TLocalModularityM() {
}

double NetworKit::TLocalModularityM::getValue(node v) {
	double inside = 0;
	double outside = 0;
	bool modified = false;
	if (community->find(v) == community->end()) {
		modified = true;
	}
	community->insert(v);
	for (node u : (*community)) {
		this->G->forNeighborsOf(u, [&](node x){
			if (community->find(x) == community->end()){
				outside ++;
			} else {
				if (u == x) {
					inside++;
				} else {
					inside = inside + 0.5;
				}
			}
		});
	}

	if (modified == true) {
		community->erase(v);
	}

	if (outside == 0) {
		return G->numberOfEdges();
	}
	return inside / outside;
}

NetworKit::TConductance::TConductance(const Graph& G, std::unordered_set<node>& community) : TQualityObjective(G, community) {
}

NetworKit::TConductance::~TConductance() {
}

double NetworKit::TConductance::getValue(node v) {

	int count = 0;
	bool modified = false;
	if (community->find(v) == community->end()) {
		modified = true;
	}
	community->insert(v);

	this->G->forNeighborsOf(v, [&](node u){
		if (community->find(u) == community->end()) {
			count++;
		}
	});

	if (modified == true) {
		community->erase(v);
	}
	if (degSum - volume - this->G->degree(v)  == 0) {
		return 0;
	}
	if (G->hasEdge(v, v)) {
		return 1 - ((double)(nBoundaryEdges + 2 * count - this->G->degree(v) + 1)/ ((double)std::min(volume + this->G->degree(v), degSum - volume - this->G->degree(v))));
	}

	return 1 - ((double)(nBoundaryEdges + 2 * count - this->G->degree(v))/ ((double)std::min(volume + this->G->degree(v), degSum - volume - this->G->degree(v))));
}
