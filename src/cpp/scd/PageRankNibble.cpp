/*
 * PageRankNibble.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include "PageRankNibble.h"
#include "ApproximatePageRank.h"
#include "../community/Conductance.h"
#include <cmath>
#include <vector>

namespace NetworKit {

PageRankNibble::PageRankNibble(Graph& g, double alpha, double epsilon): SelectiveCommunityDetector(g), alpha(alpha), epsilon(epsilon) {

}



std::set<node> PageRankNibble::bestSweepSet(const std::vector<double>& pr) {
	double entry = 0.0;
	double deg = 0.0;
	std::vector<std::pair<double, node> > sweepVec;

	// fill sweep set vector, use only non-zero entries
	G.forNodes([&](node v) {
		if (pr[v] > 0.0) {
			deg = (double) G.degree(v);
			entry = (deg > 0.0) ? (pr[v] / deg) : (0.0);
			sweepVec.push_back(std::make_pair(entry, v));
		}
	});
	count suppSize = sweepVec.size();
	TRACE("Support size: ", suppSize);


	// order vertices
	TRACE("Before sorting");
	auto comp([&](const std::pair<double, node>& a, const std::pair<double, node>& b) {
		return a.first > b.first;
	});
	std::sort(sweepVec.begin(), sweepVec.end(), comp);
	TRACE("After sorting");


	// find best sweep set w.r.t. conductance
	std::set<node> suitableCluster, bestCluster;
	Conductance conductance;
	double bestCond = std::numeric_limits<double>::max();
	double cut = 0.0;
	double volume = 0.0;

	for (index j = 0; j < suppSize; ++j) {
		// update sweep set
		node v = sweepVec[j].second;
		G.forNeighborsOf(v, [&](node neigh) {
			if (suitableCluster.count(neigh) == 0) {
				cut += G.weight(v, neigh);
			}
		});
		volume += G.degree(v);
		suitableCluster.insert(v);

		// compute conductance
		double cond = cut / volume;

		if (cond < bestCond) {
			bestCluster = suitableCluster;
		}
	}

	return bestCluster;
}


std::set<node> PageRankNibble::expandSeed(node seed) {
	DEBUG("APR(G, ", alpha, ", ", epsilon, ")");
	ApproximatePageRank apr(G, alpha, epsilon);
	std::vector<double> pr = apr.run(seed);

	std::set<node> cluster = bestSweepSet(pr);
	return cluster;
}

std::map<node, std::set<node> >  PageRankNibble::run(std::set<unsigned int>& seeds) {
    std::map<node, std::set<node> > result;
	for (auto seed : seeds) {
		auto community = expandSeed(seed);
		result[seed] = community;
	}
    return result;
}

} /* namespace NetworKit */
