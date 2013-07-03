/*
 * TGreedyCommunityExpansion.h
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#ifndef TGREEDYCOMMUNITYEXPANSION_H_
#define TGREEDYCOMMUNITYEXPANSION_H_

#include "SelectiveCommunityDetector.h"

namespace NetworKit {

template<class QualityObjective, class Acceptability, class Trimming> class TGreedyCommunityExpansion: public SelectiveCommunityDetector {

public:

	TGreedyCommunityExpansion(const Graph& G);

	virtual ~TGreedyCommunityExpansion();

	std::unordered_map<node, std::pair<std::unordered_set<node>, int64_t>> run(
			std::unordered_set<node> set);

	/**
	 * @param[in]	s	seed node
	 *
	 * @param[out]		community as a set of nodes
	 */
	std::unordered_set<node> expandSeed(node s);

protected:

	// TODO: move to postprocessing virtual double clusterClusterSimilarity (std::unordered_set<node>& community1, std::unordered_set<node>& community2);
};

template<class QualityObjective, class Acceptability, class Trimming>
inline TGreedyCommunityExpansion<QualityObjective, Acceptability, Trimming>::TGreedyCommunityExpansion(
		const Graph& G) :
		SelectiveCommunityDetector(G) {

}

template<class QualityObjective, class Acceptability, class Trimming>
inline TGreedyCommunityExpansion<QualityObjective, Acceptability, Trimming>::~TGreedyCommunityExpansion() {
}

template<class QualityObjective, class Acceptability, class Trimming>
std::unordered_map<node, std::pair<std::unordered_set<node>, int64_t>> TGreedyCommunityExpansion<
		QualityObjective, Acceptability, Trimming>::run(
		std::unordered_set<node> set) {


	std::unordered_map<node, std::pair<std::unordered_set<node>, int64_t>> communities;
	std::pair<std::unordered_set<node>, int64_t> tmp;
	Aux::Timer running;
	for (node u : set) {
		running.start();
		std::unordered_set<node> community = this->expandSeed(u);
		running.stop();
		 tmp = {community, running.elapsedMilliseconds()};
		communities.insert(std::pair<node, std::pair<std::unordered_set<node>, int64_t>>(u, tmp));
	}

	return communities;

}

template<class QualityObjective, class Acceptability, class Trimming>
inline std::unordered_set<node> TGreedyCommunityExpansion<QualityObjective,
		Acceptability, Trimming>::expandSeed(node s) {
	std::unordered_set<node> community;
	std::unordered_set<node> shell; // shell are the nodes outside of the
									// community with edges to nodes inside
	std::unordered_map<node, count> boundary;

	Acceptability acceptability(G, community, shell);
	QualityObjective objective(G, community, boundary);

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
	objective.nNodes = 1;
	objective.nBoundaryEdges = shell.size();
	objective.volume = G.degree(s);
	if (G.hasEdge(s, s)) {
		objective.nInternEdges++;
	}
	double currentObjectiveValue = objective.getValue(s)[0];
	node vMax = *(shell.begin()); // initialize vMax with a random node from the shell
	double acceptanceMax = acceptability.getValue(vMax); // maximum acceptance value

	while (expanded) {

		if (shell.empty()) {
			break;
		}
		expanded = false;
		for (node v : shell) {
			acceptanceValues.insert(
					std::pair<node, double>(v, acceptability.getValue(v)));
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
					if (objective.getValue(x)[0] - objective.getValue(vMax)[0]
							> 0.00001) {
						vMax = x;
						acceptanceMax = acc;
					} else if (objective.getValue(vMax)[0]
							- objective.getValue(x)[0] < 0.00001) {
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

			std::vector<double> tmp = objective.getValue(vMax);
			// include only nodes which lead to a strictly positive improvement
			if (tmp[0] > currentObjectiveValue) {
				expanded = true;
				currentObjectiveValue = tmp[0];
				objective.nNodes++;
				objective.nBoundaryEdges = tmp[1];
				objective.volume = objective.volume + G.degree(vMax);
				objective.nInternEdges = tmp[2];
				community.insert(vMax);
				if (G.hasEdge(vMax, vMax)) {
					boundary.insert( { vMax, G.degree(vMax) - 1 });
				} else {
					boundary.insert( { vMax, G.degree(vMax) });
				}
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

	Trimming trimm;
	trimm.run(community, G);

	return community;
}

} /* namespace NetworKit */

#endif /* TGREEDYCOMMUNITYEXPANSION_H_ */
