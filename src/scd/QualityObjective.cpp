/*
 * QualityObjective.cpp

 *
 *  Created on: 16.06.2013
 *      Author: Yassine Marrakchi
 */

#include "QualityObjective.h"

namespace NetworKit {

QualityObjective::QualityObjective(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary):nInternEdges(0), nNodes(0), nBoundaryEdges(0), volume(0) {
	this->G = &G;
	this->community = &community;
	this->boundary = &boundary;
	this->degSum = this->G->parallelSumForNodes([&](node u){
		return this->G->degree(u);
	});
}

QualityObjective::~QualityObjective() {
}

LocalModularityM::LocalModularityM(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary)
	: QualityObjective(G, community, boundary) {
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
		result.push_back(G->numberOfEdges());
		return result;
	}
	result.push_back(core/((double)boundary));
	result.push_back(boundary);
	result.push_back(core);
	return result;
}


Conductance::Conductance(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary)
		: QualityObjective(G, community, boundary) {

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
		result.push_back(0);
		result.push_back(degSum);
		return result;
	}
	if (G->hasEdge(v, v)) {
		double tmp = ((double)(nBoundaryEdges + 2 * count - this->G->degree(v) + 1)/ ((double)std::min(volume + this->G->degree(v), degSum - volume - this->G->degree(v))));
		result.push_back(1-tmp);
		result.push_back(nBoundaryEdges + 2 * count - this->G->degree(v) + 1);
		result.push_back(0);
		return result;
	}

	double tmp = (((double)(nBoundaryEdges + 2 * count - this->G->degree(v)))/ ((double)std::min(volume + this->G->degree(v), degSum - volume - this->G->degree(v))));
	result.push_back(1-tmp);
	result.push_back(nBoundaryEdges + 2 * count - this->G->degree(v));
	result.push_back(0);
	return result;
}

LocalModularityL::LocalModularityL(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary)
	: QualityObjective(G, community, boundary) {
}

LocalModularityL::~LocalModularityL() {
}

std::vector<double> LocalModularityL::getValue(node v) {

	std::vector<double> result;
	int inside = 0;
	int outside = 0;
	int bound = this->boundary->size();
	bool modified = false;
	if (community->find(v) == community->end()) {
		modified = true;
	}
	community->insert(v);
	G->forNeighborsOf(v, [&](node u) {
		if (u != v) {
			if (community->find(u) == community->end()) {
				outside ++;
			} else {
				inside++;
			}
			if (this->community->find(u) != this->community->end()) {
				if(this->boundary->find(u)->second == 1) {
					bound--;
				}
			}
		}
	});

	if (outside > 0 && modified) {
		bound++;
	}

	if (modified == true) {
		community->erase(v);
	}
	int core = inside;
	if(G->hasEdge(v, v)) {
		inside++;
	}


	if (bound == 0) {
		result.push_back(G->numberOfEdges());
		result.push_back(0);
		result.push_back(G->numberOfEdges());
		return result;
	}
	result.push_back((((double)(this->nInternEdges + inside))/(community->size() + 1))/((this->nBoundaryEdges - core + outside)/bound));
	result.push_back(this->nBoundaryEdges - core + outside);
	result.push_back(this->nInternEdges + inside);
	return result;
}

}


