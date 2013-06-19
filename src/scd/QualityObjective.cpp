/*
 * QualityObjective.cpp

 *
 *  Created on: 16.06.2013
 *      Author: Yassine Marrakchi
 */

#include "QualityObjective.h"

namespace NetworKit {

QualityObjective::QualityObjective(Graph& G, std::unordered_set<node>& community) {
	this->G = &G;
	this->community = &community;
}

QualityObjective::~QualityObjective() {
}

LocalModularityM::LocalModularityM(Graph& G, std::unordered_set<node>& community)
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


Conductance::Conductance(Graph& G, std::unordered_set<node>& community)
	: QualityObjective(G, community){
}

Conductance::~Conductance() {
}

double Conductance::getValue(node v) {
	double volume = 0;
	double boundary = 0;
	double all = 0;
	bool modified = false;
	if (community->find(v) == community->end()) {
		modified = true;
	}
	community->insert(v);

	for (node u : (*community)) {
		volume = volume + this->G->degree(u);
		this->G->forNeighborsOf(u, [&](node v){
			if (community->find(v) == community->end()) {
				boundary++;
			}
		});
	}

	G->forNodes([&](node v){
		all = all + G->degree(v);
	});

	if (modified == true) {
		community->erase(v);
	}
	if (volume == 0 || all-volume == 0)
		return 0;
	return 1 - (boundary / std::min(volume, all-volume));
}


}


