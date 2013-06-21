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
		overlapSizes[O] += 1;
	});



	/**
	 * Return the size of the union of C and D.
	 * C is a cluster in zeta, D in eta
	 */
	auto unionSize = [&](cluster C, cluster D){
		cluster inter = intersect[C][D];
		assert (inter != none); // entry must be set
		TRACE("clusters sized " << size_zeta[C] << " and " << size_eta[D] << " have overlap of " << overlapSizes[inter]);
		assert (overlapSizes[inter] >= 0);
		assert (overlapSizes[inter] <= std::max(size_zeta[C], size_eta[D]));
		return size_zeta[C] + size_eta[D] - overlapSizes[inter];
	};


	auto log_b = Aux::MissingMath::log_b; // import convenient logarithm function


	// calculate mutual information
	//		 $MI(\zeta,\eta):=\sum_{C\in\zeta}\sum_{D\in\eta}\frac{|C\cap D|}{n}\cdot\log_{2}\left(\frac{|C\cup D|\cdot n}{|C|\cdot|D|}\right)$
	double MI = 0.0; // mutual information
	for (cluster C = 0; C < zeta.upperBound(); C++) {
		for (cluster D = 0; D < eta.upperBound(); D++) {
			count sizeC = size_zeta[C];
			count sizeD = size_eta[D];
			if ((sizeC != 0) && (sizeD != 0)) {
				double factor1 = overlapSizes[intersect[C][D]] / (double) n;
				assert ((size_zeta[C] * size_eta[D]) != 0);
				TRACE("union of " << C << " and " << D << ": " << unionSize(C, D));
				double frac2 = (unionSize(C, D) * n) / (size_zeta[C] * size_eta[D]);
				double factor2 = log_b(frac2, 2);
				MI += factor1 * factor2;
			}
		}
	}

	// sanity check
	assert (! std::isnan(MI));
	assert (MI >= 0.0);
	assert (MI <= 1.0);



	// compute entropy for both clusterings
	// $H(\zeta):=-\sum_{C\in\zeta}P(C)\cdot\log_{2}(P(C))$
	double H_zeta = 0.0;
	for (cluster C = zeta.lowerBound(); C< zeta.upperBound(); ++C) {
		if (P_zeta[C] != 0) {
			H_zeta += P_zeta[C] * log_b(P_zeta[C], 2);
		} // log(0) is not defined
	}
	H_zeta = -1.0 * H_zeta;

	double H_eta = 0.0;
	for (cluster C = zeta.lowerBound(); C< zeta.upperBound(); ++C) {
		if (P_zeta[C] != 0) {
			H_eta += P_eta[C] * log_b(P_eta[C], 2);
		} // log(0) is not defined
	}
	H_eta = -1.0 * H_eta;

	assert (! std::isnan(H_zeta));
	assert (! std::isnan(H_eta));

	// entropy values range from 0 for the 1-clustering to log_2(n) for the singleton clustering
	assert (H_zeta >= 0.0);
	assert (H_eta >= 0.0);
	assert (H_zeta <= log_b(n, 2));
	assert (H_zeta <= log_b(n, 2));

	// calculate NMID
	// $NMI(\zeta,\eta):=\frac{2\cdot MI(\zeta,\eta)}{H(\zeta)+H\text{(\eta)}}$
	// $NMID(\zeta,\eta):=\begin{cases}
	//	1-NMI(\zeta,\eta)\\
	//	0 & H(\zeta)+H(\eta)=0
	//	\end{cases}$$
	double NMID;
	double NMI;
	if ((H_zeta + H_eta) == 0) {
		NMID = 0.0;
	} else {
		NMI = (2 * MI) / (H_zeta + H_eta);
		NMID = 1.0 - NMI;
	}
	// sanity check
	assert (NMID >= 0.0);
	assert (NMID <= 1.0);
	return NMID;

}

} /* namespace NetworKit */
