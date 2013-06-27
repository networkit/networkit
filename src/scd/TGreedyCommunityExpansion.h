/*
 * TGreedyCommunityExpansion.h
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#ifndef TGREEDYCOMMUNITYEXPANSION_H_
#define TGREEDYCOMMUNITYEXPANSION_H_

#include "SelectiveCommunityDetector.h"

#include "Acceptability.h"

namespace NetworKit {

template <class QualityObjective, class Acceptability, class Trimming> class TGreedyCommunityExpansion: public SelectiveCommunityDetector {

public:

	TGreedyCommunityExpansion(const Graph& G);

	virtual ~TGreedyCommunityExpansion();

	virtual std::unordered_map<node, std::unordered_set<node>> run(std::unordered_set<node> set);

	/**
	 * @param[in]	s	seed node
	 *
	 * @param[out]		community as a set of nodes
	 */
	virtual std::unordered_set<node> expandSeed(node s);

protected:

	// TODO: move to postprocessing virtual double clusterClusterSimilarity (std::unordered_set<node>& community1, std::unordered_set<node>& community2);
};



template<class QualityObjective, class Acceptability, class Trimming>
inline TGreedyCommunityExpansion<QualityObjective, Acceptability, Trimming>::TGreedyCommunityExpansion(const Graph& G) : SelectiveCommunityDetector(G) {

}

template<class QualityObjective, class Acceptability, class Trimming>
inline TGreedyCommunityExpansion<QualityObjective, Acceptability, Trimming>::~TGreedyCommunityExpansion() {
}

template<class QualityObjective, class Acceptability, class Trimming>
inline std::unordered_map<node, std::unordered_set<node> > TGreedyCommunityExpansion<QualityObjective, Acceptability, Trimming>::run(std::unordered_set<node> set) {

	std::unordered_map<node, std::unordered_set<node>> communities;

	for (node u : set) {
		std::unordered_set<node> community = this->expandSeed(u);
		communities.insert(std::pair<node, std::unordered_set<node>>(u, community));
	}

	return communities;

}

template<class QualityObjective, class Acceptability, class Trimming>
inline std::unordered_set<node> TGreedyCommunityExpansion<QualityObjective, Acceptability, Trimming>::expandSeed(node s) {
	std::unordered_set<node> community;
	std::unordered_set<node> shell; // shell are the nodes outside of the
									// community with edges to nodes inside

	std::unordered_map<node, double> acceptanceValues;
	Acceptability acceptability(G, community, shell);
	QualityObjective objective(G, community);

	double currentObjectiveValue = objective.getValue(s);
	community.insert(s);
	bool expanded = true;	// community has been expanded in the last iteration
	// all neighbors of s form the shell
	G.forNeighborsOf(s, [&](node v) {
		shell.insert(v);
	});

	if (shell.empty()) {
		return community;
	}

	objective.nBoundaryEdges = shell.size();
	objective.volume = G.degree(s);
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
					if (objective.getValue(x) - objective.getValue(vMax) > 0.00001) {
						vMax = x;
						acceptanceMax = acc;
					} else if (objective.getValue(vMax) - objective.getValue(x) < 0.00001) {
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

			if (objective.getValue(vMax) > currentObjectiveValue) {

				expanded = true;
				currentObjectiveValue = objective.getValue(vMax);
				int count = 0;
				G.forNeighborsOf(vMax, [&](node u){
					if (community.find(u) == community.end()) {
						count++;
					}
				});
				if (G.hasEdge(vMax, vMax)) {
					objective.nBoundaryEdges = objective.nBoundaryEdges + 2 * count - G.degree(vMax) + 1;
				} else {
					objective.nBoundaryEdges = objective.nBoundaryEdges + 2 * count - G.degree(vMax);
				}
				objective.volume = objective.volume + G.degree(vMax);
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

	Trimming trimm;
	trimm.run(community, G);

	return community;
}

} /* namespace NetworKit */


#endif /* TGREEDYCOMMUNITYEXPANSION_H_ */
