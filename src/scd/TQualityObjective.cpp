/*
 * TQualityObjective.cpp
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#include "TQualityObjective.h"

namespace NetworKit {

TQualityObjective::TQualityObjective(const Graph& G,
		std::unordered_set<node>& community,
		std::unordered_map<node, count>& boundary) :
		nInternEdges(0), nNodes(0), nBoundaryEdges(0), volume(0) {
	this->G = &G;
	this->community = &community;
	this->boundary = &boundary;
	this->degSum = this->G->parallelSumForNodes([&](node u) {
		return this->G->degree(u);
	});
}

TQualityObjective::~TQualityObjective() {
}

std::vector<double> TQualityObjective::getValue(node v) {
	throw std::runtime_error(
			"abstract method called - must be redefined in subclass");
}

TLocalModularityM::TLocalModularityM(const Graph& G,
		std::unordered_set<node>& community,
		std::unordered_map<node, count>& boundary) :
		TQualityObjective(G, community, boundary) {
}

TLocalModularityM::~TLocalModularityM() {
}

std::vector<double> TLocalModularityM::getValue(node v) {

	std::vector<double> result;
	count inside = 0;
	count outside = 0;
	bool modified = false;
	if (community->find(v) == community->end()) {
		modified = true;
	}
	if (!modified) {
		result.push_back(
				((double) (this->volume - this->nBoundaryEdges))
						/ this->nBoundaryEdges);
		result.push_back(0);
		result.push_back(0);
		return result;
	}
	community->insert(v);
	this->G->forNeighborsOf(v, [&](node x) {
		if (community->find(x) == community->end()) {
			outside ++;
		} else {
			if (x != v) {
				inside++;
			}
		}
	});

	community->erase(v);
	count boundary = this->nBoundaryEdges - inside + outside;
	if (G->hasEdge(v, v)) {
		inside++;
	}
	count core = this->nInternEdges + inside;
	if (boundary == 0) {
		result.push_back(G->numberOfEdges());
		result.push_back(0);
		result.push_back(G->numberOfEdges());
		return result;
	}
	result.push_back(core / ((double) boundary));
	result.push_back(boundary);
	result.push_back(core);
	return result;
}

TConductance::TConductance(const Graph& G, std::unordered_set<node>& community,
		std::unordered_map<node, count>& boundary) :
		TQualityObjective(G, community, boundary) {

}

TConductance::~TConductance() {
}

std::vector<double> TConductance::getValue(node v) {

	count count = 0;
	bool modified = false;
	std::vector<double> result;
	if (community->find(v) == community->end()) {
		modified = true;
	}

	community->insert(v);

	this->G->forNeighborsOf(v, [&](node u) {
		if (community->find(u) == community->end()) {
			count++;
		}
	});

	if (modified == true) {
		community->erase(v);
	}
	if (degSum - volume - this->G->degree(v) == 0) {
		result.push_back(0);
		result.push_back(0);
		result.push_back(degSum);
		return result;
	}
	if (G->hasEdge(v, v)) {
		double tmp = ((double) (nBoundaryEdges + 2 * count - this->G->degree(v)
				+ 1)
				/ ((double) std::min(volume + this->G->degree(v),
						degSum - volume - this->G->degree(v))));
		result.push_back(1 - tmp);
		result.push_back(nBoundaryEdges + 2 * count - this->G->degree(v) + 1);
		result.push_back(0);
		return result;
	}

	double tmp = (((double) (nBoundaryEdges + 2 * count - this->G->degree(v)))
			/ ((double) std::min(volume + this->G->degree(v),
					degSum - volume - this->G->degree(v))));
	result.push_back(1 - tmp);
	result.push_back(nBoundaryEdges + 2 * count - this->G->degree(v));
	result.push_back(0);
	return result;
}

TLocalModularityL::TLocalModularityL(const Graph& G,
		std::unordered_set<node>& community,
		std::unordered_map<node, count>& boundary) :
		TQualityObjective(G, community, boundary) {
}

TLocalModularityL::~TLocalModularityL() {
}

std::vector<double> TLocalModularityL::getValue(node v) {
	std::vector<double> result;
	count inside = 0; 	// number of ingoing edges from v to community
	count outside = 0;	// number of outgoing edges from v to commmunity
	count bound = this->boundary->size();
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
	count core = inside;
	if (G->hasEdge(v, v)) {
		inside++;
	}

	if (bound == 0) {
		result.push_back(G->numberOfEdges());
		result.push_back(0);
		result.push_back(G->numberOfEdges());
		return result;
	}
	result.push_back(
			(((double) (this->nInternEdges + inside)) / (community->size() + 1))
					/ ((this->nBoundaryEdges - core + outside) / bound));
	result.push_back(this->nBoundaryEdges - core + outside);
	result.push_back(this->nInternEdges + inside);
	return result;
}

} /* namespace NetworKit */

