/*
 * NormalizedMutualInformation.cpp
 *
 *  Created on: 30.04.2013
 *      Author: cls
 */

#include "NMIDistance.h"

#include "DynamicNMIDistance.h"
#include "../auxiliary/MissingMath.h"
#include "../auxiliary/NumericTools.h"
#include "../auxiliary/Log.h"
#include "PartitionIntersection.h"

namespace NetworKit {


double NMIDistance::getDissimilarity(const Graph& G, const Partition& zeta, const Partition& eta) {

	count n = G.numberOfNodes();

	DEBUG("zeta=" , zeta.getVector());
	DEBUG("eta=" , eta.getVector());


	std::vector<count> size_zeta(zeta.upperBound());
	std::vector<count> size_eta(eta.upperBound());

	// precompute sizes for each cluster

	G.forNodes([&](node u){
		index C = zeta[u];
		index D = eta[u];
		assert (C != none);
		assert (D != none);
		size_zeta[C] += 1;
		size_eta[D] += 1;
	});

	DEBUG("size_zeta=" , size_zeta);
	DEBUG("size_eta=" , size_eta);

	// precompute cluster probabilities
	std::vector<double> P_zeta(zeta.upperBound(), 0.0);
	std::vector<double> P_eta(eta.upperBound(), 0.0);

	for (index C = zeta.lowerBound(); C < zeta.upperBound(); ++C) {
		P_zeta[C] = size_zeta[C] / (double) n;
	}

	for (index C = eta.lowerBound(); C < eta.upperBound(); ++C) {
		P_eta[C] = size_eta[C] / (double) n;
	}


	Partition overlap = PartitionIntersection().calculate(zeta, eta);
	DEBUG("overlap=", overlap.getVector());

	std::vector<std::vector<index> > intersect(zeta.upperBound(), std::vector<index>(eta.upperBound(), none)); // intersect[C][D] returns the overlap cluster

	std::vector<count> overlapSizes(overlap.upperBound(), 0); // overlapSizes[O] returns the size of the overlap cluster

	G.forNodes([&](node u){
		index O = overlap[u];
		index C = zeta[u];
		index D = eta[u];
		intersect[C][D] = O; // writes are redundant - but this may be efficient

	});

	G.forNodes([&](node u){
		index O = overlap[u];
		overlapSizes[O] += 1;
	});

	DEBUG("overlapSizes=", overlapSizes);



	/**
	 * Return the size of the union of C and D.
	 * C is a cluster in zeta, D in eta
	 */
#if (LOG_LEVEL != LOG_LEVEL_INFO)
	auto unionSize = [&](index C, index D){
		index O = intersect[C][D];
		assert (O != none); // entry must be set
		count sizeC = size_zeta[C];
		count sizeD = size_eta[D];
		count sizeO = overlapSizes[O];

		TRACE("clusters sized " , sizeC , " and " , sizeD , " have overlap of " , sizeO);
		assert (sizeO >= 0);
		assert (sizeO <= std::max(sizeC, sizeD));
		return sizeC + sizeD - sizeO;
	};
#endif


	auto log_b = Aux::MissingMath::log_b; // import convenient logarithm function


	// calculate mutual information
	//		 $MI(\zeta,\eta):=\sum_{C\in\zeta}\sum_{D\in\eta}\frac{|C\cap D|}{n}\cdot\log_{2}\left(\frac{|C\cap D|\cdot n}{|C|\cdot|D|}\right)$
	double MI = 0.0; // mutual information
	for (index C = 0; C < zeta.upperBound(); C++) {
		for (index D = 0; D < eta.upperBound(); D++) {
			count sizeC = size_zeta[C];
			count sizeD = size_eta[D];
			if ((sizeC != 0) && (sizeD != 0)) { // cluster ids may not correspond to a real cluster
				// the two clusters may or may not intersect
				index O = intersect[C][D];
				if (O == none) { // clusters do not intersect
					TRACE("clusters do not intersect: " , C , ", " , D);
				} else {
					count sizeO = overlapSizes[O];
					double factor1 =  sizeO / (double) n;
					assert ((size_zeta[C] * size_eta[D]) != 0);
					TRACE("overlap of " , C , " and " , D , " has size: " , sizeO);
					TRACE("union of " , C , " and " , D , " has size: " , unionSize(C, D));
					double frac2 = (sizeO * n) / (double) (size_zeta[C] * size_eta[D]);
					assert (frac2 != 0);
					double factor2 = log_b(frac2, 2);
					TRACE("frac2 = " , frac2 , ", factor1 = " , factor1 , ", factor2 = " , factor2);
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
