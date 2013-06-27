/*
 * DynamicNMIDistance.cpp
 *
 *  Created on: Jun 26, 2013
 *      Author: Henning
 */

#include "DynamicNMIDistance.h"

namespace NetworKit {

DynamicNMIDistance::DynamicNMIDistance() {
	// TODO Auto-generated constructor stub

}

DynamicNMIDistance::~DynamicNMIDistance() {
	// TODO Auto-generated destructor stub
}

double DynamicNMIDistance::getDissimilarity(Graph& newGraph,
		Clustering& oldClustering, Clustering& newClustering) {

	auto isInBoth = [&](node u, const Clustering& oldClustering, const Clustering& newClustering) {
		return ((newClustering[u] != none) && (u < oldClustering.numberOfEntries()));
	};

	count n = newGraph.numberOfNodes();

	DEBUG("oldClustering=" << Aux::vectorToString(oldClustering.getVector()));
	DEBUG("newClustering=" << Aux::vectorToString(newClustering.getVector()));


	std::vector<count> size_old(oldClustering.upperBound());
	std::vector<count> size_new(newClustering.upperBound());

	// precompute sizes for each cluster
	newGraph.forNodes([&](node u){
		if (isInBoth(u, oldClustering, newClustering)) {
			cluster C = oldClustering[u];
			cluster D = newClustering[u];
			size_old[C]++;
			size_new[D]++;
		}
	});

	DEBUG("size_old=" << Aux::vectorToString(size_old));
	DEBUG("size_new=" << Aux::vectorToString(size_new));

	// precompute cluster probabilities
	std::vector<double> P_old(oldClustering.upperBound(), 0.0);
	std::vector<double> P_new(newClustering.upperBound(), 0.0);

	double n_double = (double) n;
	for (cluster C = oldClustering.lowerBound(); C < oldClustering.upperBound(); ++C) {
		P_old[C] = size_old[C] / n_double;
	}
	for (cluster C = newClustering.lowerBound(); C < newClustering.upperBound(); ++C) {
		P_new[C] = size_new[C] / n_double;
	}


	HashingOverlapper hashing;
	std::vector<Clustering> clusterings;
	clusterings.push_back(oldClustering);
	clusterings.push_back(newClustering);

	Clustering overlap = hashing.run(newGraph, clusterings);
	DEBUG("overlap=" << Aux::vectorToString(overlap.getVector()));

	std::vector<std::vector<cluster> > intersect(oldClustering.upperBound(),
			std::vector<cluster>(newClustering.upperBound(), none)); // intersect[C][D] returns the overlap cluster

	std::vector<count> overlapSizes(overlap.upperBound(), 0); // overlapSizes[O] returns the size of the overlap cluster

	newGraph.forNodes([&](node u){
		if (isInBoth(u, oldClustering, newClustering)) {
			cluster O = overlap[u];
			cluster C = oldClustering[u];
			cluster D = newClustering[u];
			intersect[C][D] = O; // writes are redundant - but this may be efficient
		}
	});

	newGraph.forNodes([&](node u){
		if (isInBoth(u, oldClustering, newClustering)) {
			cluster O = overlap[u];
			overlapSizes[O]++;
		}
	});

	DEBUG("overlapSizes=" << Aux::vectorToString(overlapSizes));



	/**
	 * Return the size of the union of C and D.
	 * C is a cluster in zeta, D in eta
	 */
	auto unionSize = [&](cluster C, cluster D){
		cluster O = intersect[C][D];
		assert (O != none); // entry must be set
		count sizeC = size_old[C];
		count sizeD = size_new[D];
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
	for (cluster C = 0; C < oldClustering.upperBound(); C++) {
		for (cluster D = 0; D < newClustering.upperBound(); D++) {
			count sizeC = size_old[C];
			count sizeD = size_new[D];
			if ((sizeC != 0) && (sizeD != 0)) { // cluster ids may not correspond to a real cluster
				// the two clusters may or may not intersect
				cluster O = intersect[C][D];
				if (O == none) { // clusters do not intersect
					TRACE("clusters do not intersect: " << C << ", " << D);
				} else {
					count sizeO = overlapSizes[O];
					double factor1 =  sizeO / (double) n;
					assert ((size_old[C] * size_new[D]) != 0);
					TRACE("union of " << C << " and " << D << " has size: " << unionSize(C, D));
					TRACE("overlap of " << C << " and " << D << " has size: " << sizeO);
					double frac2 = (sizeO * n) / (double) (size_old[C] * size_new[D]);
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
	double H_old = entropy(oldClustering, n, P_old);
	double H_new = entropy(newClustering, n, P_new);

	// calculate NMID
	// $NMI(\zeta,\eta):=\frac{2\cdot MI(\zeta,\eta)}{H(\zeta)+H\text{(\eta)}}$
	// $NMID(\zeta,\eta):=\begin{cases}
	//	1-NMI(\zeta,\eta)\\
	//	0 & H(\zeta)+H(\eta)=0
	//	\end{cases}$$
	double NMID = 0.0;
	double NMI = 0.0;
	double H_sum = H_old + H_new;
	combineValues(H_sum, MI, NMI, NMID);
	sanityCheck(NMI, NMID);

	return NMID;
}

void DynamicNMIDistance::combineValues(double H_sum, double MI, double& NMI, double& NMID) const {
	if (Aux::NumericTools::equal(H_sum, 0.0)) {
		NMID = 0.0;
	} else {
		NMI = (2.0 * MI) / H_sum;
		NMID = 1.0 - NMI;
	}
}

double DynamicNMIDistance::entropy(const Clustering& clustering, count n, std::vector<double> probs) {
	auto log_b = Aux::MissingMath::log_b; // import convenient logarithm function

	// $H(\zeta):=-\sum_{C\in\zeta}P(C)\cdot\log_{2}(P(C))$
	double H = 0.0;
	for (cluster C = clustering.lowerBound(); C < clustering.upperBound(); ++C) {
		if (probs[C] != 0) {
			H += probs[C] * log_b(probs[C], 2);
		} // log(0) is not defined
	}
	H = -1.0 * H;

	assert (! std::isnan(H));

	// entropy values range from 0 for the 1-clustering to log_2(n) for the singleton clustering
	assert (H >= 0.0);
	assert (H <= log_b(n, 2));

	return H;
}

void DynamicNMIDistance::sanityCheck(double& NMI, double& NMID) const {
	assert (Aux::NumericTools::ge(NMI, 0.0));
	assert (Aux::NumericTools::le(NMI, 1.0));

	// if NMID is close to 0 because of numerical error
	DEBUG("sanity check, NMID: " << NMID);
	if (Aux::NumericTools::equal(NMID, 0.0)) {
		NMID = 0.0;
	}
	TRACE("sanity check, NMID after: " << NMID);

	assert (Aux::NumericTools::ge(NMID, 0.0));
	assert (Aux::NumericTools::le(NMID, 1.0));
}


} /* namespace NetworKit */
