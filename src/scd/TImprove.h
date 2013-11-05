/*
 * TImprove.h
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#ifndef TIMPROVE_H_
#define TIMPROVE_H_

#include "SelectiveCommunityDetector.h"

namespace NetworKit {

template<class QualityObjective> class TImprove: public SelectiveCommunityDetector {

public:

	TImprove(const Graph& G);

	virtual ~TImprove();

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

template<class QualityObjective>
inline TImprove<QualityObjective>::TImprove(
		const Graph& G) :
		SelectiveCommunityDetector(G) {
}

template<class QualityObjective>
inline TImprove<QualityObjective>::~TImprove() {
}

template<class QualityObjective>
std::unordered_map<node, std::pair<std::unordered_set<node>, int64_t>> TImprove<
		QualityObjective>::run(
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

template<class QualityObjective>
inline std::unordered_set<node> TImprove<QualityObjective>::expandSeed(node s) {

	std::unordered_set<node> community;
	std::unordered_set<node> shell; // shell are the nodes outside of the
									// community with edges to nodes inside
	std::unordered_map<node, count> boundary;
	std::unordered_set<node> candidates;
	node current;
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

	objective.nNodes = 1;
	objective.nBoundaryEdges = shell.size();
	objective.volume = G.degree(s);
	if (G.hasEdge(s, s)) {
		objective.nInternEdges++;
	}
	double currentObjectiveValue = objective.getValue(s)[0];
	while (expanded) {

		if (shell.empty()) {
			break;
		}
		candidates.clear();
		for (node u : shell) {
			candidates.insert(u);
		}
		expanded = false;

		while (!expanded && !candidates.empty()) {
			current = *candidates.cbegin();
			std::vector<double> tmp = objective.getValue(current);
			if (tmp[0] > currentObjectiveValue) {
				expanded = true;
				currentObjectiveValue = tmp[0];
				objective.nNodes++;
				objective.nBoundaryEdges = tmp[1];
				objective.volume = objective.volume + G.degree(current);
				objective.nInternEdges = tmp[2];
				community.insert(current);
				if (G.hasEdge(current, current)) {
					boundary.insert( { current, G.degree(current) - 1 });
				} else {
					boundary.insert( { current, G.degree(current) });
				}
				shell.erase(current);
				candidates.erase(current);
				G.forNeighborsOf(current, [&](node v) {
					if (community.find(v) == community.end()) {
						shell.insert(v);
						candidates.insert(v);
					} else {
						if(boundary.find(v)->second == 1) {
							boundary.erase(boundary.find(v));
						} else {
							boundary.find(v)->second--;
						}
					}
				});
			} else {
				candidates.erase(current);
			}
		} // end while shell.size() != 0
	} // end while expanded

	return community;
}

} /* namespace NetworKit */

#endif /* TIMPROVE_H_ */
