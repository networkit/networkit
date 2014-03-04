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

PageRankNibble::PageRankNibble(Graph& g): G(g) {

}

PageRankNibble::~PageRankNibble() {

}


std::set<node> PageRankNibble::suitableSweepSet(const std::vector<double>& pr) {
	count upperNodeId = G.upperNodeIdBound();
	count suppSize = this->supportSize(pr);
	DEBUG("Support size: ", suppSize);
//	double totalEdgeWeight = G.totalEdgeWeight();
//	double volThrsh = 4.0 / 3.0 * totalEdgeWeight;
	double deg = 0.0;

	std::vector<std::pair<double, node> > sweepVec(upperNodeId);
	G.forNodes([&](node v) {
		deg = (double) G.degree(v);
		sweepVec[v].first = (deg > 0.0) ? (pr[v] / deg) : (0.0);
		sweepVec[v].second = v;
	});


	auto print([&]() {
		for (auto entry: sweepVec) {
			std::cout << entry.first << " ";
		}
		std::cout << std::endl;
	});


	// order vertices, use only supportSize many afterwards
	DEBUG("Before sorting");
	std::sort(sweepVec.begin(), sweepVec.end());
	// reverse
	std::reverse(sweepVec.begin(), sweepVec.end());
	// TODO: directly sort in descending order instead
	DEBUG("After sorting");


	// find best sweep set w.r.t. conductance
	std::set<node> suitableCluster, bestCluster;
	count volume = 0;
	Partition partition(upperNodeId);
	partition.allToOnePartition();
	Conductance conductance;
	double bestCond = std::numeric_limits<double>::max();
	node first = sweepVec[0].second;
	partition.toSingleton(first);
	index id = partition[first];

	for (index j = 0; j < suppSize; ++j) {
		// update sweep set
		node v = sweepVec[j].second;
		partition.moveToSubset(id, v);
		volume += G.degree(v);
		suitableCluster.insert(v);

		// compute conductance
		double cond = conductance.getQuality(partition, G); // TODO: accelerate

		if (cond < bestCond) {
			bestCluster = suitableCluster;
		}
	}

	return bestCluster;
}

count PageRankNibble::supportSize(const std::vector<double>& vec) const {
	count size = 0;

	// count non-zero vector entries, tolerate numerical errors (TODO: check if disadvantageous)
	for (auto entry: vec) {
		if (entry > 0.0) {
			++size;
		}
	}

	return size;
}

std::set<node> PageRankNibble::run(node seed, double alpha, double epsilon) {
	DEBUG("APR(G, ", alpha, ", ", epsilon, ")");
	ApproximatePageRank apr(G, alpha, epsilon);
	std::vector<double> pr = apr.run(seed);

	std::set<node> cluster = suitableSweepSet(pr);
	return cluster;
}

} /* namespace NetworKit */
