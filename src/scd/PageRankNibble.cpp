/*
 * PageRankNibble.cpp
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#include "PageRankNibble.h"
#include "ApproximatePageRank.h"
#include "../clustering/Conductance.h"
#include <cmath>
#include <vector>

namespace NetworKit {

PageRankNibble::PageRankNibble(Graph& g): G(g) {

}

PageRankNibble::~PageRankNibble() {

}

std::set<node> PageRankNibble::suitableSweepSet(const std::vector<double>& pr, double phi, unsigned int b, unsigned int B) {
	count n = G.numberOfNodes();
	count suppSize = this->supportSize(pr);
	std::vector<std::pair<double, node> > sweepVec;

	G.forNodes([&](node v) {
		sweepVec[v].first = pr[v] / G.degree(v);
		sweepVec[v].second = v;
	});

	// order vertices, use only supportSize many afterwards
	std::sort(sweepVec.begin(), sweepVec.end());


	// check conditions
	auto isSetSuitable([&](double cond, double vol) {
		bool result = cond < phi;

		if (result) {
			// volume condition
			result = vol > (1 << (b-1));
			result = result && (vol < 4.0 * G.totalEdgeWeight() / 3.0);
		}

		if (result) {
			// probability change condition (similar to spectral gap)
			index pos1 = 1 << b;
			index pos2 = 1 << (b-1);
			result = pr[pos1] - pr[pos2] > 1.0 / (48.0 * B);
		}

		return result;
	});


	// for each possible set: check conditions
	std::set<node> suitableCluster;
	count volume = 0;
	Partition partition(n, 0);
	Conductance conductance;

	for (index j = 0; j < suppSize; ++j) {
		// update sweep set
		node v = sweepVec[j].second;
		partition[v] = 1;
		volume += G.degree(v);
		suitableCluster.insert(v);

		// compute conductance
		double cond = conductance.getQuality(partition, G);

		// check conditions
		if (isSetSuitable(cond, volume)) {
			return suitableCluster;
		}
	}

	// return empty set if not fulfilled for any set
	suitableCluster.clear();
	return suitableCluster;
}

count PageRankNibble::supportSize(const std::vector<double>& vec) const {
	count size = 0;

	// count non-zero vector entries, tolerate numerical errors (TODO: check if disadvantageous)
	for (auto entry: vec) {
		if (! Aux::NumericTools::equal(entry, 0.0)) {
			++size;
		}
	}

	return size;
}

std::set<node> PageRankNibble::run(node seed, double phi, unsigned int b) {
	count m = G.numberOfEdges();
	count B = ceil(log2(m));
	double alpha = phi * phi / (225.0 * log(100.0 * sqrt(m)));
	double epsilon = pow(2, -b) / (48.0 * B);

	DEBUG("APR(G, ", alpha, ", ", epsilon);
	ApproximatePageRank apr(G, alpha, epsilon);
	std::vector<double> pr = apr.run(seed);

	std::set<node> cluster = suitableSweepSet(pr, phi, b, B);
	return cluster;
}

} /* namespace NetworKit */
