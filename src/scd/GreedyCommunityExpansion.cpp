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
	INFO("running GreedyCommunityExpansion from node " << s);

	std::unordered_set<node> community;
	std::unordered_set<node> shell; // shell are the nodes outside of the
	                                // community with edges to nodes inside

	std::map<node, double> acceptanceValues;

	// TODO: make these selectable later
	GreedyCommunityExpansion::DummySimilarity acceptability(G, community, shell);
	GreedyCommunityExpansion::Conductance conductance(G, community);

	double currentObjectiveValue = conductance.getValue(s);


	community.insert(s);
	TRACE("adding start node " << s);
	bool expanded = true;		// community has been expanded in the last iteration

	// all neighbors of s form the shell
	G.forNeighborsOf(s, [&](node v) {
		shell.insert(v);
	});



	if (shell.empty()) {
		return community;
	}

	node vMax = *(shell.begin()); // initialize vMax with a random node from the shell
	double acceptanceMax = acceptability.getValue(vMax);	// maximum acceptance value

	while (expanded) {
		if (shell.empty()) {
			break;
		}
		expanded = false;
		for (node v : shell) {
			acceptanceValues.insert(std::pair<node,double> (v, acceptability.getValue(v)));
		}
		while (!acceptanceValues.empty()) {
			acceptanceMax = 0;
			for (auto it = acceptanceValues.begin(); it != acceptanceValues.end(); ++it ) {
				node x = it->first;
				double acc = it->second;
				if (it->second > acceptanceMax) {
					vMax = x;
					acceptanceMax = acc;
				} else if (it ->second == acceptanceMax) {
					if (conductance.getValue(vMax) < conductance.getValue(x)) {
						vMax = x;
						acceptanceMax = acc;
					} else if (conductance.getValue(vMax) == conductance.getValue(x)) {
						// first tie-breaking by degree
						if (G.degree(x) > G.degree(vMax)) {
							vMax = x;
							acceptanceMax = acc;
						} else if (G.degree(x) == G.degree(vMax)) {
							// last tie-breaking by id
							if (x < vMax) {
								vMax = it->first;
								acceptanceMax = it->second;
							}
						}
					}
				}

			}

			// include only nodes which lead to a strictly positive improvement
			if (conductance.getValue(vMax) > currentObjectiveValue) {
				expanded = true;
				currentObjectiveValue = conductance.getValue(vMax);
				community.insert(vMax);
				TRACE("adding node " << vMax);
				shell.erase(vMax);
				G.forNeighborsOf(vMax, [&](node v){
					if (community.find(v) == community.end()) {
						shell.insert(v);
					}
				});
				acceptanceValues.clear();
			} else {
				auto it = acceptanceValues.find(vMax);
				node x = it->first;
				// node with highest acceptability is discarded from the map
				acceptanceValues.erase(x);

			}
		} // end while acceptanceValues.size() != 0

	} // end while expanded

	DummyTrimming trimm;
	trimm.run(community, G);

	return community;
}


GreedyCommunityExpansion::QualityObjective::QualityObjective(
	Graph& G, std::unordered_set<node>& community) {
	this->G = &G;
	this->community = &community;
}

GreedyCommunityExpansion::LocalModularityM::LocalModularityM(Graph& G, std::unordered_set<node>& community)
	: QualityObjective(G, community){
}

GreedyCommunityExpansion::LocalModularityM::~LocalModularityM() {
}

double GreedyCommunityExpansion::LocalModularityM::getValue(node v) {

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

GreedyCommunityExpansion::Acceptability::Acceptability(
	Graph& G, std::unordered_set<node>& community, std::unordered_set<node>& shell) {
	this->G = &G;
	this->community = &community;
	this->shell = &shell;
}


GreedyCommunityExpansion::NodeClusterSimilarity::NodeClusterSimilarity(
	Graph& G, std::unordered_set<node>& community, std::unordered_set<node>& shell): Acceptability(G, community, shell) {
}

GreedyCommunityExpansion::NodeClusterSimilarity::~NodeClusterSimilarity() {
}

double GreedyCommunityExpansion::NodeClusterSimilarity::getValue(node v) {

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

GreedyCommunityExpansion::QualityObjective::~QualityObjective() {
}

GreedyCommunityExpansion::Acceptability::~Acceptability() {
}

GreedyCommunityExpansion::Conductance::Conductance(Graph& G, std::unordered_set<node>& community)
	: QualityObjective(G, community){
}

GreedyCommunityExpansion::Conductance::~Conductance() {
}

double GreedyCommunityExpansion::Conductance::getValue(node v) {
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

GreedyCommunityExpansion::DummySimilarity::DummySimilarity(
		Graph& G, std::unordered_set<node>& community, std::unordered_set<node>& shell): Acceptability(G, community, shell) {
}

GreedyCommunityExpansion::DummySimilarity::~DummySimilarity() {
}

double GreedyCommunityExpansion::DummySimilarity::getValue(node v) {
	return 0.5;
}

std::map<node, std::unordered_set<node> > GreedyCommunityExpansion::seedSetExpansion(
		Graph& G, std::vector<node> set) {

	std::map<node, std::unordered_set<node>> communities;
	GreedyCommunityExpansion GCE;
	for (node u : set) {
		communities.insert(std::pair<node, std::unordered_set<node>> (u, GCE.run(G, u)));
	}
	bool modified = true;

	return communities;
}

double GreedyCommunityExpansion::clusterClusterSimilarity(
		std::unordered_set<node>& community1,
		std::unordered_set<node>& community2) {
	double intersection = 0;
	for (node u : community1) {
		if (community2.find(u) != community2.end()) {
			intersection++;
		}
	}
	return intersection / (community1.size() + community2.size() - intersection);
}

} /* namespace NetworKit */


