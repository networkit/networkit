/*
 * GreedyCommunityExpansion.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include "GreedyCommunityExpansion.h"

namespace NetworKit {

GreedyCommunityExpansion::GreedyCommunityExpansion(const Graph& G,
		Acceptability& similarity, QualityObjective& objective,
		CommunityTrimming& trimming) :
		SelectiveCommunityDetector(G) {
	this->objective = &objective;
	this->similarity = &similarity;
	this->trimming = &trimming;
}

GreedyCommunityExpansion::~GreedyCommunityExpansion() {
	// TODO Auto-generated destructor stub
}

std::unordered_set<node> GreedyCommunityExpansion::expandSeed(node s) {

	std::unordered_set<node> community;
	std::unordered_set<node> shell; // shell are the nodes outside of the
									// community with edges to nodes inside
	std::unordered_map<node, count> boundary;
	this->similarity->community = &community;
	this->similarity->shell = &shell;
	this->objective->community = &community;
	this->objective->boundary = &boundary;

	std::unordered_map<node, double> acceptanceValues;
	community.insert(s);
	if (G.hasEdge(s, s)) {
		boundary.insert( { s, G.degree(s) - 1 });
	} else {
		boundary.insert( { s, G.degree(s) });
	}
	bool expanded = true;	// community has been expanded in the last iteration

	// all neighbors of s form the shell
	G.forNeighborsOf(s, [&](node v) {
		shell.insert(v);
	});

	if (shell.empty()) {
		return community;
	}
	objective->nNodes = 1;
	objective->nBoundaryEdges = shell.size();
	objective->volume = G.degree(s);
	if (G.hasEdge(s, s)) {
		objective->nInternEdges++;
	}
	double currentObjectiveValue = objective->getValue(s)[0];

	node vMax = *(shell.begin()); // initialize vMax with a random node from the shell
	double acceptanceMax = similarity->getValue(vMax); // maximum acceptance value

	while (expanded) {
		if (shell.empty()) {
			break;
		}
		expanded = false;
		for (node v : shell) {
			acceptanceValues.insert(
					std::pair<node, double>(v, similarity->getValue(v)));
		}

		while (!acceptanceValues.empty()) {
			acceptanceMax = 0;
			for (auto it = acceptanceValues.begin();
					it != acceptanceValues.end(); ++it) {

				node x = it->first;
				double acc = it->second;
				if (it->second - acceptanceMax > 0.00001) {
					vMax = x;
					acceptanceMax = acc;
				} else if (it->second - acceptanceMax < 0.00001) {
					if (objective->getValue(x)[0] - objective->getValue(vMax)[0]
							> 0.00001) {
						vMax = x;
						acceptanceMax = acc;
					} else if (objective->getValue(vMax)[0]
							- objective->getValue(x)[0] < 0.00001) {
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

			std::vector<double> tmp = objective->getValue(vMax);

			// include only nodes which lead to a strictly positive improvement
			if (tmp[0] > currentObjectiveValue) {

				expanded = true;
				currentObjectiveValue = tmp[0];
				objective->nNodes++;
				objective->nBoundaryEdges = tmp[1];
				objective->volume = objective->volume + G.degree(vMax);
				objective->nInternEdges = tmp[2];
				community.insert(vMax);
				boundary.insert( { vMax, G.degree(vMax) });
				shell.erase(vMax);
				G.forNeighborsOf(vMax, [&](node v) {
					if (community.find(v) == community.end()) {
						shell.insert(v);
					} else {
						if(boundary.find(v)->second == 1) {
							boundary.erase(boundary.find(v));
						} else {
							boundary.find(v)->second--;
						}
					}
				});
				acceptanceValues.clear();
			} else {
				// node with highest acceptability is discarded from the map
				acceptanceValues.erase(acceptanceValues.find(vMax)->first);
			}
		} // end while acceptanceValues.size() != 0
	} // end while expanded

	trimming->run(community, G);
	return community;
}

std::unordered_map<node, std::unordered_set<node>> GreedyCommunityExpansion::run(
		std::unordered_set<node> set) {

	std::unordered_map<node, std::unordered_set<node>> communities;

	for (node u : set) {
		std::unordered_set<node> community = this->expandSeed(u);
		communities.insert(
				std::pair<node, std::unordered_set<node>>(u, community));
	}

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

