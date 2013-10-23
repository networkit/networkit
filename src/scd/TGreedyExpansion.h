/* TGreedyExpansion.h
 *
 *  Created on: 14.09.2013
 *      Author: cls, Yassine Marrakchi
 */

#ifndef TGREEDYEXPANSION_H_
#define TGREEDYEXPANSION_H_

#include "SelectiveCommunityDetector.h"

namespace NetworKit {

template<class QualityObjective, class Similarity> class TGreedyExpansion: public SelectiveCommunityDetector {

public:

	TGreedyExpansion(const Graph& G);

	virtual ~TGreedyExpansion();

	std::unordered_map<node, std::pair<std::unordered_set<node>, int64_t>> run(
			std::unordered_set<node> set);

	/**
	 * @param[in]	s	seed node
	 *
	 * @param[out]		community as a set of nodes
	 */
	std::unordered_set<node> expandSeed(node s);

protected:

};

template<class QualityObjective, class Similarity>
inline TGreedyExpansion<QualityObjective, Similarity>::TGreedyExpansion(
		const Graph& G) :
		SelectiveCommunityDetector(G) {

}

template<class QualityObjective, class Similarity>
inline TGreedyExpansion<QualityObjective, Similarity>::~TGreedyExpansion() {
}

template<class QualityObjective, class Similarity>
std::unordered_map<node, std::pair<std::unordered_set<node>, int64_t>> TGreedyExpansion<QualityObjective,
	Similarity>::run(std::unordered_set<node> set) {

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

template<class QualityObjective, class Similarity>
inline std::unordered_set<node> TGreedyExpansion<QualityObjective, Similarity>::expandSeed(node s) {

	std::unordered_set<node> community;
	std::unordered_set<node> shell; // shell are the nodes outside of the
									// community with edges to nodes inside
	std::unordered_map<node, count> boundary;
	std::unordered_map<node, double> similarityValues;

	Similarity similarity(G);
	QualityObjective objective(G, community, boundary);

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
	double product = 0;
	double temporary0 = 0;
	double temporary1 = 0;
	double temporary2 = 0;
	node v = s;


	while (expanded) {
		if (shell.empty()) {
				break;
		}
		expanded = false;

		G.forNeighborsOf(v, [&](node x) {
			if (community.find(v) == community.end()) {
				shell.insert(v);
			} else {
				if(boundary.find(v)->second == 1) {
					boundary.erase(boundary.find(v));
				} else {
					boundary.find(v)->second--;
				}
			}
			if (similarityValues.find(x) != similarityValues.end()) {
				temporary2 = similarity.getValue(x, s);
				if (temporary2 > similarityValues.find(x)->second) {
					similarityValues.find(x)->second = temporary2;
				}
			} else {
				similarityValues.insert(std::pair<node, double>(x, similarity.getValue(x, s)));
			}
		});
		for (node u : shell) {
			temporary0 = objective.getValue(u)[0];
			if (temporary0 > currentObjectiveValue) {
				temporary1 = similarityValues.find(u)->second;
				if (temporary0 * temporary1 > product) {
					v = u;
					currentObjectiveValue = temporary0;
					product = temporary0 * temporary1;
				}
			}
		}
		if (community.find(v) == community.end()) {
			std::vector<double> tmp = objective.getValue(v);
			expanded = true;
			currentObjectiveValue = tmp[0];
			objective.nNodes++;
			objective.nBoundaryEdges = tmp[1];
			objective.volume = objective.volume + G.degree(v);
			objective.nInternEdges = tmp[2];
			community.insert(v);
			if (G.hasEdge(v, v)) {
				boundary.insert( { v, G.degree(v) - 1 });
			} else {
				boundary.insert( { v, G.degree(v) });
			}
			shell.erase(v);
		}
	}
	return community;
}

} /* namespace NetworKit */

#endif /* TGREEDYEXPANSION_H_ */
