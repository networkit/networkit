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

	DEBUG("zeta=" << Aux::vectorToString(zeta.getVector()));
	DEBUG("eta=" << Aux::vectorToString(eta.getVector()));


	std::vector<count> size_zeta(zeta.upperBound());
	std::vector<count> size_eta(eta.upperBound());

	// precompute sizes for each cluster

	G.forNodes([&](node u){
		cluster C = zeta[u];
		cluster D = eta[u];
		assert (C != none);
		assert (D != none);
		size_zeta[C] += 1;
		size_eta[D] += 1;
	});

	DEBUG("size_zeta=" << Aux::vectorToString(size_zeta));
	DEBUG("size_eta=" << Aux::vectorToString(size_eta));

	// precompute cluster probabilities
	std::vector<double> P_zeta(zeta.upperBound(), 0.0);
	std::vector<double> P_eta(eta.upperBound(), 0.0);

	for (cluster C = zeta.lowerBound(); C < zeta.upperBound(); ++C) {
		P_zeta[C] = size_zeta[C] / (double) n;
	}

	for (cluster C = eta.lowerBound(); C < eta.upperBound(); ++C) {
		P_eta[C] = size_eta[C] / (double) n;
	}



	// RegionGrowingOverlapper hashing;
	HashingOverlapper hashing;
	std::vector<Clustering> clusterings;
	clusterings.push_back(zeta);
	clusterings.push_back(eta);

	Clustering overlap = hashing.run(G, clusterings);
	DEBUG("overlap=" << Aux::vectorToString(overlap.getVector()));

	std::vector<std::vector<cluster> > intersect(zeta.upperBound(), std::vector<cluster>(eta.upperBound(), none)); // intersect[C][D] returns the overlap cluster

	std::vector<count> overlapSizes(overlap.upperBound(), 0); // overlapSizes[O] returns the size of the overlap cluster

	G.forNodes([&](node u){
		cluster O = overlap[u];
		cluster C = zeta[u];
		cluster D = eta[u];
		intersect[C][D] = O; // writes are redundant - but this may be efficient

	});

	G.forNodes([&](node u){
		cluster O = overlap[u];
		overlapSizes[O] += 1;
	});

	DEBUG("overlapSizes=" << Aux::vectorToString(overlapSizes));



	/**
	 * Return the size of the union of C and D.
	 * C is a cluster in zeta, D in eta
	 */
	auto unionSize = [&](cluster C, cluster D){
		cluster O = intersect[C][D];
		assert (O != none); // entry must be set
		count sizeC = size_zeta[C];
		count sizeD = size_eta[D];
		count sizeO = overlapSizes[O];

		TRACE("clusters sized " << sizeC << " and " << sizeD << " have overlap of " << sizeO);
		assert (sizeO >= 0);
		assert (sizeO <= std::max(sizeC, sizeD));
		return sizeC + sizeD - sizeO;
	};


	auto log_b = Aux::MissingMath::log_b; // import convenient logarithm function


	// calculate mutual information
	//		 $MI(\zeta,\eta):=\sum_{C\in\zeta}\sum_{D\in\eta}\frac{|C\cap D|}{n}\cdot\log_{2}\left(\frac{|C\cup D|\cdot n}{|C|\cdot|D|}\right)$
	double MI = 0.0; // mutual information
	for (cluster C = 0; C < zeta.upperBound(); C++) {
		for (cluster D = 0; D < eta.upperBound(); D++) {
			count sizeC = size_zeta[C];
			count sizeD = size_eta[D];
			if ((sizeC != 0) && (sizeD != 0)) { // cluster ids may not correspond to a real cluster
				// the two clusters may or may not intersect
				cluster O = intersect[C][D];
				if (O == none) { // clusters do not intersect
					TRACE("clusters do not intersect: " << C << ", " << D);
				} else {
					count sizeO = overlapSizes[O];
					double factor1 =  sizeO / (double) n;
					assert ((size_zeta[C] * size_eta[D]) != 0);
					TRACE("union of " << C << " and " << D << " has size: " << unionSize(C, D));
					TRACE("overlap of " << C << " and " << D << " has size: " << sizeO);
					double frac2 = (sizeO * n) / (double) (size_zeta[C] * size_eta[D]);
					assert (frac2 != 0);
					double factor2 = log_b(frac2, 2);
					TRACE("frac2 = " << frac2 << ", factor1 = " << factor1 << ", factor2 = " << factor2);
					MI += factor1 * factor2;
				}
			}
		}
	}

	// sanity check
	assert (! std::isnan(MI));
	assert (MI >= 0.0);

	// compute entropy for both clusterings
	DynamicNMIDistance dynNMID;
	double H_zeta = dynNMID.entropy(zeta, n, P_zeta);
	double H_eta = dynNMID.entropy(eta, n, P_eta);

	assert (! std::isnan(H_zeta));
	assert (! std::isnan(H_eta));

	double H_sum = H_zeta + H_eta;
	double NMI = 0.0;
	double NMID = 0.0;
	dynNMID.combineValues(H_sum, MI, NMI, NMID);
	dynNMID.sanityCheck(NMI, NMID);

	return NMID;
}

} /* namespace NetworKit */
