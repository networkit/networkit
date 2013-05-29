/*
 * NormalizedMutualInformation.cpp
 *
 *  Created on: 30.04.2013
 *      Author: cls
 */

#include "NMIDistance.h"

namespace NetworKit {

NMIDistance::NMIDistance() {
	// TODO Auto-generated constructor stub

}

NMIDistance::~NMIDistance() {
	// TODO Auto-generated destructor stub
}

double NMIDistance::getDissimilarity(Graph& G, Clustering& zeta, Clustering& eta) {

	count n = G.numberOfNodes();

	std::vector<count> zetaSizes(zeta.upperBound());
	std::vector<count> etaSizes(eta.upperBound());

	// calculate sizes for each cluster
	zeta.forEntries([&](node u, cluster C){
		zetaSizes[C] += 1;
	});

	eta.forEntries([&](node u, cluster C){
		etaSizes[C] += 1;
	});

	HashingOverlapper hashing;
	std::vector<Clustering> clusterings;
	clusterings.push_back(zeta);
	clusterings.push_back(eta);

	Clustering overlap = hashing.run(G, clusterings);



	cluster kZeta = zeta.upperBound();
	cluster kEta = eta.upperBound();

	std::vector<std::vector<cluster> > intersect(kZeta, std::vector<cluster>(kEta, none)); // intersection[C][D] returns O, the overlap cluster

	std::vector<count> overlapSizes(overlap.upperBound(), 0); // overlap[O] returns the size of the overlap cluster

	G.forNodes([&](node u){
		cluster O = overlap[u];
		cluster C = zeta[u];
		cluster D = eta[u];
		intersect[C][D] = O; // writes are redundant - but this may be efficient
		overlapSizes[O] += 1;
	});


	/**
	 * log_b(x)
	 */
	auto logb = [](double x, double b){
		return log(x) / log(b);
	};


	/**
	 * Return the size of the union of C and D.
	 * C is a cluster in zeta, D in eta
	 */
	auto unionSize = [&](cluster C, cluster D){
		return zetaSizes[C] + etaSizes[D] - overlapSizes[intersect[C][D]];
	};

	double mi = 0.0; // mutual information
	for (cluster C = 0; C < kZeta; C++) {
		for (cluster D = 0; D < kEta; D++) {
			double factor1 = overlapSizes[intersect[C][D]] / (double) n;
			double factor2 = logb((unionSize(C, D) * n) / (zetaSizes[C] * etaSizes[D]), 2);
			mi += factor1 * factor2;
		}
	}

	// sanity check
	assert (mi >= 0.0);
	assert (mi <= 1.0);

	auto size = [&](cluster C, Clustering& zeta) {
		return 0; // TODO: implement
	};

	auto P = [&](cluster C, Clustering& zeta) {
		return size(C, zeta) / (double) n;
	};

	auto H = [&](Clustering& zeta) {
		// TODO: implement
	};

	auto MI = [&](Clustering& zeta, Clustering& eta) {

	};
}

} /* namespace NetworKit */
