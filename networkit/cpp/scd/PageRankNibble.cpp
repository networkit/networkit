/*
 * PageRankNibble.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include "PageRankNibble.h"
#include "ApproximatePageRank.h"
#include "../community/Conductance.h"
#include "../auxiliary/Parallel.h"
#include <cmath>
#include <vector>

namespace NetworKit {

PageRankNibble::PageRankNibble(Graph& g, double alpha, double epsilon): SelectiveCommunityDetector(g), alpha(alpha), epsilon(epsilon) {
	assert(!g.isWeighted());
}

PageRankNibble::~PageRankNibble() {

}


std::set<node> PageRankNibble::bestSweepSet(std::vector<std::pair<node, double>>& pr) {
	TRACE("Support size: ", pr.size());


	// order vertices
	TRACE("Before sorting");
	auto comp([&](const std::pair<node, double>& a, const std::pair<node, double>& b) {
		return (a.second / G.degree(a.first)) > (b.second / G.degree(b.first));
	});
	Aux::Parallel::sort(pr.begin(), pr.end(), comp);
	TRACE("After sorting");

	for (std::vector<std::pair<node, double>>::iterator it = pr.begin(); it != pr.end(); it++) {
		TRACE("(", it->first, ", ", it->second, ")");
	}


	// find best sweep set w.r.t. conductance
	double bestCond = std::numeric_limits<double>::max();
	double cut = 0.0;
	double volume = 0.0;
	index bestSweepSetIndex = 0;
	std::unordered_map<node, bool> withinSweepSet;
	std::vector<node> currentSweepSet;

	for (std::vector<std::pair<node, double>>::iterator it = pr.begin(); it != pr.end(); it++) {
		// update sweep set
		node v = it->first;
		G.forNeighborsOf(v, [&](node neigh) {
			if (withinSweepSet.find(neigh) == withinSweepSet.end()) {
				cut++;
			} else {
				cut--;
			}
		});
		volume += G.volume(v);
		currentSweepSet.push_back(v);
		withinSweepSet[v] = true;

		// compute conductance
		double cond = cut / fmin(volume, 2 * G.numberOfEdges() - volume);

		std::stringstream debug;

		debug << "Current vertex: " << v << "; Current sweep set conductance: " << cond << std::endl;
		debug << "Current cut weight: " << cut << "; Current volume: " << volume << std::endl;
		debug << "Total graph volume: " << 2 * G.numberOfEdges() << std::endl;

		TRACE(debug.str());

		if (cond < bestCond) {
			bestCond = cond;
			bestSweepSetIndex = currentSweepSet.size();
		}
	}

	std::set<node> bestSweepSet;

	for (index j = 0; j < bestSweepSetIndex; j++) {
		bestSweepSet.insert(currentSweepSet[j]);
	}

	return bestSweepSet;
}


std::set<node> PageRankNibble::expandSeed(node seed) {
	DEBUG("APR(G, ", alpha, ", ", epsilon, ")");
	ApproximatePageRank apr(G, alpha, epsilon);
	std::vector<std::pair<node, double>> pr = apr.run(seed);
	std::set<node> s = bestSweepSet(pr);
	return s;
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
