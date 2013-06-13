/*
 * GreedyCommunityExpansion.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include "GreedyCommunityExpansion.h"

namespace NetworKit {

GreedyCommunityExpansion::GreedyCommunityExpansion() {
	// TODO Auto-generated constructor stub

}

GreedyCommunityExpansion::~GreedyCommunityExpansion() {
	// TODO Auto-generated destructor stub
}

std::unordered_set<node> GreedyCommunityExpansion::run(Graph& G, node s) {

	std::unordered_set<node> community;
	std::unordered_set<node> shell; // shell are the nodes outside of the
	                                // community with edges to nodes inside

	std::map<node,double,double> acceptanceValues;

	GreedyCommunityExpansion::NodeClusterSimilarity acceptability(G, &community, &shell);
	GreedyCommunityExpansion::Conductance conductance(G, &community);

	double currentObjectiveValue = conductance.getValue(s);
	node vMax;
	double acceptanceMax = 0;

	community.insert(s);

	G.forNeighborsOf(s, [&](node v) {
		shell.insert(v);
	});


	do {
		for (node v : shell) {
		}







		if(conductance.getValue(vMax) <= currentObjectiveValue){
			community.insert(vMax);
			shell.erase(vMax);
			G.forNeighborsOf(vMax, [&](node v){
				if (community.find(v) == community.end()) {
					shell.insert(v);
				}
			});
		}
	} while (true);


	// TODO: optional trimming phase according to node fitness
	return community;
}


GreedyCommunityExpansion::QualityObjective::QualityObjective(
	Graph& G, std::unordered_set<node>* community): G(G), community(community) {
}

GreedyCommunityExpansion::LocalModularityM::LocalModularityM(Graph& G, std::unordered_set<node>* community)
	: QualityObjective(G, community){
}

GreedyCommunityExpansion::LocalModularityM::~LocalModularityM() {
}

double GreedyCommunityExpansion::LocalModularityM::getValue(node v) {

	double inside = 0;
	double outside = 0;
	for (auto it = community->begin(); it != community->end(); ++it) {
		this->G.forNeighborsOf(*it, [&](node v){
			if (community->find(v) == community->end()){
				outside ++;
			} else {
				if (*it == v) {
					inside++;
				} else {
					inside = inside + 0.5;
				}
			}
		});
	}

	return inside / outside;
}

GreedyCommunityExpansion::Acceptability::Acceptability(
	Graph& G, std::unordered_set<node>* community, std::unordered_set<node>* shell): G(G), community(community), shell(shell) {
}


GreedyCommunityExpansion::NodeClusterSimilarity::NodeClusterSimilarity(
	Graph& G, std::unordered_set<node>* community, std::unordered_set<node>* shell): Acceptability(G, community, shell) {
}

GreedyCommunityExpansion::NodeClusterSimilarity::~NodeClusterSimilarity() {
}

double GreedyCommunityExpansion::NodeClusterSimilarity::getValue(node v) {

	int intersection = 0;
	this->G.forNeighborsOf(v, [&](node u){
	//	if ((*(this->community)).find(u) < (*(this->community)).end()||(*(this->shell)).find(u) < (*(this->shell)).end()) { // fixme: there is a fault hier
//			intersection++;
	//	}
	});
	if (G.hasEdge(v, v)) {
		return intersection / (G.degree(v) + (*community).size() + (*shell).size() - intersection);
	} else {
		return (intersection + 1) / (G.degree(v) + (*community).size() + (*shell).size() - intersection);
	}
	return 0;
}

GreedyCommunityExpansion::QualityObjective::~QualityObjective() {
}

GreedyCommunityExpansion::Acceptability::~Acceptability() {
}

GreedyCommunityExpansion::Conductance::Conductance(Graph& G, std::unordered_set<node>* community)
	: QualityObjective(G, community){
}

GreedyCommunityExpansion::Conductance::~Conductance() {
}

double GreedyCommunityExpansion::Conductance::getValue(node v) {
	double volume = 0;
	double boundary = 0;
	double all = 0;
	community->insert(v);

	for (auto it = community->begin(); it != community->end(); ++it) {
		volume = volume + this->G.degree(*it);
		this->G.forNeighborsOf(*it, [&](node v){
			if (community->find(v) == community->end()) boundary++;
		});
	}

	G.forNodes([&](node v){
		all = all + G.degree(v);
	});
	community->erase(v);
	if (volume == 0 || all-volume == 0)
		return 1;
	return boundary / std::min(volume, all-volume);
}

} /* namespace NetworKit */
