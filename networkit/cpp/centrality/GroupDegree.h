/*
 * GroupDegree.h
 *
 *  Created on: 20.04.2018
 *      Author: Eugenio Angriman
 */

#ifndef GROUPDEGREE_H_
#define GROUPDEGREE_H_

#include "../auxiliary/BucketPQ.h"
#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup centrality
 */
class GroupDegree : public Algorithm {

public:
	/**
	 * Finds the group with the highest group degree centrality according to the
	 * definition proposed in 'The centrality of groups and classes' by Everett et
	 * al. (The Journal of mathematical sociology, 1999). This is a submodular but
	 * non monotone function so the algorithm can find a solution that is at least
	 * 1/2 of the optimum. Worst-case running time is quadratic, but usually
	 * faster in real-world networks.
	 * The 'countGroupNodes' option also count the nodes inside the group in the
	 * score, this make the group degree monotone and submodular and the algorithm
	 * is guaranteed to return a (1 - 1/e)-approximation of the optimal solution.
	 *
	 * @param G A graph.
	 * @param k Size of the group of nodes
	 * @param countGroupNodes if nodes inside the group should be counted in the
	 * centrality score.
	 */
	GroupDegree(const Graph &G, count k = 1, bool countGroupNodes = true);

	/**
	 * Computes the group with maximum degree centrality of the graph passed in
	 * the constructor.
	 */
	void run() override;

	/**
	 * Returns the group with maximum degree centrality.
	 */
	std::vector<node> groupMaxDegree();

	/**
	 * Returns the score of the group with maximum degree centrality (i.e. the
	 * number of nodes outside the group that can be reached in one hop from at
	 * least one node in the group).
	 */
	count getScore();

	/**
	 * Returns the score of the given group.
	 */
	count scoreOfGroup(const std::vector<node> &group) const;

protected:
	Graph G;
	const count k;
	const bool countGroupNodes;
	count n;
	std::vector<node> group;
	std::vector<int64_t> gain;
	std::vector<bool> reachable;
	std::vector<bool> affected;
	std::vector<bool> inGroup;
	Aux::BucketPQ queue;
	count groupScore;

	void init();
	void updateQueue();
	void updateGroup();
	void computeScore();
	void checkHasRun();
	void checkGroup(const std::vector<node> &group) const;
};

inline std::vector<node> GroupDegree::groupMaxDegree() {
	checkHasRun();
	return group;
}

inline count GroupDegree::getScore() {
	checkHasRun();
	return groupScore;
}

inline void GroupDegree::computeScore() {
	groupScore = std::count(reachable.begin(), reachable.end(), true);

	if (!countGroupNodes) {
		groupScore -= k;
	}
}

inline void GroupDegree::checkHasRun() {
	if (!hasRun) {
		throw std::runtime_error("Run method has not been called.");
	}
}

inline void GroupDegree::checkGroup(const std::vector<node> &group) const {
	std::vector<node> sortedV(group);
	std::sort(sortedV.begin(), sortedV.end());
	node u;
	auto checkNode = [&](node u) {
		if (!G.hasNode(u)) {
			std::stringstream err;
			err << "Error: node" << u << " is not in the graph.";
			throw std::runtime_error(err.str());
		}
	};
	for (count i = 0; i < sortedV.size() - 1; ++i) {
		u = sortedV[i];
		checkNode(u);
		if (u == sortedV[i + 1]) {
			throw std::runtime_error("Error: the set contains duplicate elements.");
		}
	}
	checkNode(sortedV.back());
}

inline count GroupDegree::scoreOfGroup(const std::vector<node> &group) const {
	checkGroup(group);
	std::vector<bool> touched(n, false);
	for (count i = 0; i < group.size(); ++i) {
		touched[group[i]] = true;
	}

	for (count i = 0; i < group.size(); ++i) {
		node u = group[i];
		G.forNeighborsOf(u, [&](node v) {
			if (!touched[v]) {
				touched[v] = true;
			}
		});
	}

	count result = std::count(touched.begin(), touched.end(), true);
	if (!countGroupNodes) {
		result -= k;
	}

	return result;
}

} // namespace NetworKit

#endif
