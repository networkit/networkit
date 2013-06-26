/*
 * QualityObjective.cpp

 *
 *  Created on: 16.06.2013
 *      Author: Yassine Marrakchi
 */

#include "QualityObjective.h"

namespace NetworKit {

QualityObjective::QualityObjective(const Graph& G, std::unordered_set<node>& community) {
	this->G = &G;
	this->community = &community;
}

QualityObjective::~QualityObjective() {
}

LocalModularityM::LocalModularityM(const Graph& G, std::unordered_set<node>& community)
	: QualityObjective(G, community){
}

LocalModularityM::~LocalModularityM() {
}

double LocalModularityM::getValue(node v) {

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


Conductance::Conductance(const Graph& G, std::unordered_set<node>& community) : QualityObjective(G, community), degSum(0), nBoundaryEdges(0), volume(0) {
	// TODO: precomputation of degree sum should not happen more than once for a graph
	this->degSum = this->G->parallelSumForNodes([&](node u){
		return this->G->degree(u);
	});
}

Conductance::~Conductance() {
}

double Conductance::getValue(node v) {

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

LocalModularityL::LocalModularityL(const Graph& G, std::unordered_set<node>& community)
	: QualityObjective(G, community){
}

LocalModularityL::~LocalModularityL() {
}

double LocalModularityL::getValue(node v) {
	double inside = 0;
	double outside = 0;
	std::unordered_set<node> boundary;
	bool modified = false;
	if (community->find(v) == community->end()) {
		modified = true;
	}
	community->insert(v);

	for (node u : *community) {
		this->G->forNeighborsOf(u, [&](node x){
			if (community->find(x) == community->end()){
				outside ++;
				if (boundary.find(u) == boundary.end()) {
					boundary.insert(u);
				}
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
	return (inside / community->size()) / (outside / boundary.size());
}

}


