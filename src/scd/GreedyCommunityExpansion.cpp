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

std::unordered_set<node> GreedyCommunityExpansion::expandSeed(Graph& G, node s) {

	std::unordered_set<node> community;
	std::unordered_set<node> shell; // shell are the nodes outside of the
									// community with edges to nodes inside

	std::unordered_map<node, double> acceptanceValues;

	// TODO: make these selectable later
	DummySimilarity acceptability(G, community, shell);
	Conductance conductance(G, community);

	double currentObjectiveValue = conductance.getValue(s);
	int tmp = 0;

	community.insert(s);
	bool expanded = true;	// community has been expanded in the last iteration

	// all neighbors of s form the shell
	G.forNeighborsOf(s, [&](node v) {
		shell.insert(v);
	});

	if (shell.empty()) {
		return community;
	}

	node vMax = *(shell.begin()); // initialize vMax with a random node from the shell
	double acceptanceMax = acceptability.getValue(vMax);// maximum acceptance value

	while (expanded) {
		if (shell.empty()) {
			break;
		}
		expanded = false;
		for (node v : shell) {
			acceptanceValues.insert(std::pair<node, double>(v, acceptability.getValue(v)));
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
					tmp = conductance.getValue(vMax);
					if (conductance.getValue(x) - tmp > 0.00001) {

						vMax = x;
						acceptanceMax = acc;
						tmp = conductance.getValue(x);
					} else if (conductance.getValue(vMax) - conductance.getValue(x) < 0.00001) {
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
				int count = 0;
				G.forNeighborsOf(vMax, [&](node u){
					if (community.find(u) == community.end()) {
						count++;
					}
				});
				if (G.hasEdge(vMax, vMax)) {
					conductance.nBoundaryEdges = conductance.nBoundaryEdges + 2 * count - G.degree(vMax) + 1;
				} else {
					conductance.nBoundaryEdges = conductance.nBoundaryEdges + 2 * count - G.degree(vMax);
				}
				conductance.volume = conductance.volume + G.degree(vMax);
				community.insert(vMax);
				shell.erase(vMax);
				G.forNeighborsOf(vMax, [&](node v) {
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

std::unordered_map<node, std::unordered_set<node>> GreedyCommunityExpansion::run(
		Graph& G, std::unordered_set<node> set) {

	std::unordered_map<node, std::unordered_set<node>> communities;
	GreedyCommunityExpansion GCE;
	for (node u : set) {

		std::unordered_set<node> community = GCE.expandSeed(G, u);
		std::cout<<"-----------------------------"<<std::endl;


		communities.insert(std::pair<node, std::unordered_set<node>>(u, community));
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

