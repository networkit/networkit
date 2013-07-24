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

	std::vector<double> result(3, 0);
	int inside = 0;
	int outside = 0;
	bool modified = false;
	if (community->find(v) == community->end()) {
		modified = true;
	}
	if (!modified) {
		result[0] = ((double)(this->volume - this->nBoundaryEdges)) / this->nBoundaryEdges;
		result[1] = 0;
		result[2] = 0;
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
		result[0] = G->numberOfEdges();
		result[1] = 0;
		result[2] = G->numberOfEdges();
		return result;
	}
	result[0] = ((double)core)/((double)boundary);
	result[1] = boundary;
	result[2] = core;
	return result;
}


ConductanceDistance::ConductanceDistance(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary)
		: QualityObjective(G, community, boundary) {

}

ConductanceDistance::~ConductanceDistance() {
}

std::vector<double> ConductanceDistance::getValue(node v) {

	int count = 0;
	bool modified = false;
	std::vector<double> result (3, 0.0);
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
		result[0] = 0;
		result[1] = 0;
		result[2] = degSum;
		return result;
	}
	if (G->hasEdge(v, v)) {
		double tmp = ((double)(nBoundaryEdges + 2 * count - this->G->degree(v) + 1)/ ((double)std::min(volume + this->G->degree(v), degSum - volume - this->G->degree(v))));
		result[0] = 1-tmp;
		result[1] = nBoundaryEdges + 2 * count - this->G->degree(v) + 1;
		result[2] = 0;
		return result;
	}

	double tmp = (((double)(nBoundaryEdges + 2 * count - this->G->degree(v)))/ ((double)std::min(volume + this->G->degree(v), degSum - volume - this->G->degree(v))));
	result[0] = 1-tmp;
	result[1] = nBoundaryEdges + 2 * count - this->G->degree(v);
	result[2] = 0;
	return result;
}

LocalModularityL::LocalModularityL(const Graph& G, std::unordered_set<node>& community, std::unordered_map<node,count>& boundary)
	: QualityObjective(G, community, boundary) {
}

LocalModularityL::~LocalModularityL() {
}

std::vector<double> LocalModularityL::getValue(node v) {

	std::vector<double> result(3, 0);
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
		result[0] = G->numberOfEdges();
		result[1] = 0;
		result[2] = G->numberOfEdges();
		return result;
	}

	result[0] = (((double)(this->nInternEdges + inside))/(community->size() + 1))/((this->nBoundaryEdges - core + outside)/(double)bound);
	result[1] = this->nBoundaryEdges - core + outside;
	result[2] = this->nInternEdges + inside;
	return result;
}

}


