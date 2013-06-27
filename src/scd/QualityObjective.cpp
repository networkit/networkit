/*
 * QualityObjective.cpp

 *
 *  Created on: 16.06.2013
 *      Author: Yassine Marrakchi
 */

#include "QualityObjective.h"

namespace NetworKit {

QualityObjective::QualityObjective(const Graph& G, std::unordered_set<node>& community):nInternEdges(0), nNodes(0), nBoundaryEdges(0), volume(0), nBoundaryNodes(0) {
	this->G = &G;
	this->community = &community;
	this->degSum = this->G->parallelSumForNodes([&](node u){
		return this->G->degree(u);
	});
}

QualityObjective::~QualityObjective() {
}

LocalModularityM::LocalModularityM(const Graph& G, std::unordered_set<node>& community)
	: QualityObjective(G, community) {
}

LocalModularityM::~LocalModularityM() {
}

std::vector<double> LocalModularityM::getValue(node v) {

	std::vector<double> result;
	int inside = 0;
	int outside = 0;
	bool modified = false;
	if (community->find(v) == community->end()) {
		modified = true;
	}
	if (!modified) {
		result.push_back(((double)(this->volume - this->nBoundaryEdges)) / this->nBoundaryEdges);
		result.push_back(0);
		result.push_back(0);
		result.push_back(0);
		return result;
	}
	community->insert(v);
	this->G->forNeighborsOf(v, [&](node x){
		if (community->find(x) == community->end()){
			outside ++;
		} else {
			if (x != v) {
				inside++;
			}
		}
	});

	community->erase(v);
	int	boundary = this->nBoundaryEdges - inside + outside;
	if(G->hasEdge(v, v)) {
		inside++;
	}
	int core = this->nInternEdges + inside;
	if (boundary == 0) {
		result.push_back(G->numberOfEdges());
		result.push_back(0);
		result.push_back(0);
		result.push_back(0);
		return result;
	}
	result.push_back(core/((double)boundary));
	result.push_back(boundary);
	result.push_back(0);
	result.push_back(0);
	result.push_back(core);
	return result;

}


Conductance::Conductance(const Graph& G, std::unordered_set<node>& community) : QualityObjective(G, community) {

}

Conductance::~Conductance() {
}

std::vector<double> Conductance::getValue(node v) {

	int count = 0;
	bool modified = false;
	std::vector<double> result;
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
		result.push_back(0);
		result.push_back(degSum);
		result.push_back(0);
		result.push_back(0);
		return result;
	}
	if (G->hasEdge(v, v)) {
		double tmp = ((double)(nBoundaryEdges + 2 * count - this->G->degree(v) + 1)/ ((double)std::min(volume + this->G->degree(v), degSum - volume - this->G->degree(v))));
		result.push_back(1-tmp);
		result.push_back(nBoundaryEdges + 2 * count - this->G->degree(v) + 1);
		result.push_back(0);
		result.push_back(0);
		return result;
	}

	double tmp = ((double)(nBoundaryEdges + 2 * count - this->G->degree(v))/ ((double)std::min(volume + this->G->degree(v), degSum - volume - this->G->degree(v))));
	result.push_back(1-tmp);
	result.push_back(nBoundaryEdges + 2 * count - this->G->degree(v));
	result.push_back(0);
	result.push_back(0);
	return result;
}

LocalModularityL::LocalModularityL(const Graph& G, std::unordered_set<node>& community)
	: QualityObjective(G, community) {
}

LocalModularityL::~LocalModularityL() {
}

std::vector<double> LocalModularityL::getValue(node v) {
	std::vector<double> result;
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
	return result;

	//return (inside / community->size()) / (outside / boundary.size());
}

}


