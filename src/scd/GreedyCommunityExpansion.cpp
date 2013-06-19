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

	std::map<node, double> acceptanceValues;

	// TODO: make these selectable later
	GreedyCommunityExpansion::DummySimilarity acceptability(G, community, shell);
	GreedyCommunityExpansion::Conductance conductance(G, community);

	double currentObjectiveValue = conductance.getValue(s);


	community.insert(s);
	bool expanded = true;		// community has been expanded in the last iteration

	// all neighbors of s form the shell
	G.forNeighborsOf(s, [&](node v) {
		shell.insert(v);
	});

	if (shell.empty()) {
		INFO("shell of {s}�is empty because s is an isolated node");
		// if shell is empty, just return the singleton community of s
		return community;
	}

	node vMax = *(shell.begin()); // initialize vMax with a random node from the shell
	double acceptanceMax = acceptability.getValue(vMax);	// maximum acceptance value


	while(expanded) {
		if (shell.empty()) {
			INFO("there are no more nodes in the shell - breaking iteration now");
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
					if (conductance.getValue(vMax) > conductance.getValue(x)) {
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


	INFO("community size at end of run: " << community.size());
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
	for (node u : (*community)) {
		this->G->forNeighborsOf(u, [&](node v){
			if (community->find(v) == community->end()){
				outside ++;
			} else {
				if (u == v) {
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
	std::cout << this->shell->size() << std::endl;
	std::cout << this->community->size() << std::endl;


	int intersection = 0;
	this->G->forNeighborsOf(v, [&](node u) {

		if (this->community->find(u) != this->community->end()||this->shell->find(u) != this->shell->end()) {
			intersection++;
		}
	});
	if (G->hasEdge(v, v)) {
		return intersection / (G->degree(v) + (*community).size() + (*shell).size() - intersection);
	} else {


		return (intersection + 1) / (G->degree(v) + (*community).size() + (*shell).size() - intersection);
	}
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
	community->erase(v);
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

} /* namespace NetworKit */


