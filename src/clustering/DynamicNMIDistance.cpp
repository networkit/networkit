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


bool DynamicNMIDistance::isInBoth(node u, const Clustering& oldClustering, const Clustering& newClustering) {
	return ((newClustering[u] != none) &&
			(u < oldClustering.numberOfEntries()) &&
			(oldClustering[u] != none));
}


/**
 * TODO: Needs cleanup!
 * Formula follows Dhillon, Guan, Kulis: A Unified View of Kernel k-means, ...
 */
double DynamicNMIDistance::getDissimilarity(Graph& newGraph,
		Clustering& oldClustering, Clustering& newClustering) {

	auto log_b = Aux::MissingMath::log_b; // import convenient logarithm function

//	count n = newGraph.numberOfNodes();

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


	// confusion matrix
	std::vector<std::vector<count> > confMatrix = this->confusionMatrix(newGraph, oldClustering, newClustering);

	auto numOverlap = [&](Matrix& confMatrix) {
		count num = 0;
		for (Matrix::iterator iter = confMatrix.begin(); iter != confMatrix.end(); ++iter) {
			for (index i = 0; i < iter->size(); ++i) {
				num += (*iter)[i];
			}
		}
		return num;
	};

	 count totalOverlap = numOverlap(confMatrix);
	 double numDouble = (double) totalOverlap;

	/**
	 * compute \sum_{i=1}^k n_i^{(D)} with k = number of clusters in first
	 */
	auto compAggregate1 = [&](const Matrix& confMatrix, const Clustering& first, cluster D) {
		count sum = 0;
		count upperId = first.upperBound();

		for (index i = 0; i < upperId; ++i) {
			sum += confMatrix[i][D];
		}

		return sum;
	};

	/**
	 * compute \sum_{i=1}^c n_C^{(i)} with c = number of clusters in second
	 */
	auto compAggregate2 = [&](const Matrix& confMatrix, cluster C, const Clustering& second) {
		count sum = 0;
		count upperId = second.upperBound();

		for (index i = 0; i < upperId; ++i) {
			sum += confMatrix[C][i];
		}

		return sum;
	};

	double MI = 0.0; // mutual information
	for (cluster C = 0; C < oldClustering.upperBound(); C++) {
		for (cluster D = 0; D < newClustering.upperBound(); D++) {
			count currOverlap = confMatrix[C][D];
			if (currOverlap > 0) {
				double factor1 = (double) currOverlap / (double) numDouble;
				double nominator = (double) (currOverlap * numDouble);
				double aggregate1 = (double) size_old[C]; //  compAggregate1(confMatrix, oldClustering, D);
				double aggregate2 = (double) size_new[D]; // compAggregate2(confMatrix, C, newClustering);
				double denom = aggregate1 * aggregate2;
				DEBUG("frac: " << nominator << " / " << denom << " = " << nominator / denom);
				double factor2 = log_b(nominator / denom, 2);
				MI += factor1 * factor2;

				DEBUG("contribution of " << C << " and " << D << ": " << factor1 << " * " << factor2 << " = " << factor1 * factor2);
			}
		}
	}


	// precompute cluster probabilities
	std::vector<double> P_old(oldClustering.upperBound(), 0.0);
	std::vector<double> P_new(newClustering.upperBound(), 0.0);

	for (cluster C = oldClustering.lowerBound(); C < oldClustering.upperBound(); ++C) {
		P_old[C] = ((double) size_old[C]) / numDouble;
	}
	for (cluster C = newClustering.lowerBound(); C < newClustering.upperBound(); ++C) {
		P_new[C] = ((double) size_new[C]) / numDouble;
	}


#if 0


	HashingOverlapper overlapper;
	std::vector<Clustering> clusterings;
	clusterings.push_back(oldClustering);
	clusterings.push_back(newClustering);

	Clustering overlap = overlapper.run(newGraph, clusterings);
	DEBUG("overlap=" << Aux::vectorToString(overlap.getVector()));

	std::vector<std::vector<cluster> > intersect(oldClustering.upperBound(),
			std::vector<cluster>(newClustering.upperBound(), none)); // intersect[C][D] returns the overlap cluster

	std::vector<count> overlapSizes(overlap.upperBound(), 0); // overlapSizes[O] returns the size of the overlap cluster
	TRACE("upperBound = " << overlap.upperBound());

	newGraph.forNodes([&](node u){
		if (isInBoth(u, oldClustering, newClustering)) {
			cluster O = overlap[u];
			cluster C = oldClustering[u];
			cluster D = newClustering[u];
			intersect[C][D] = O; // writes are redundant - but this may be efficient
		}
	});

	count totalOverlap = 0;
	newGraph.forNodes([&](node u){
		if (isInBoth(u, oldClustering, newClustering)) {
			cluster O = overlap[u];
			TRACE("O = " << O);
			overlapSizes[O]++;
			++totalOverlap;
		}
	});
	const double totalOverlapDouble = (double) totalOverlap;

	DEBUG("overlapSizes=" << Aux::vectorToString(overlapSizes));




	double MI = 0.0;



	// calculate mutual information
	//	correct (cap): $MI(\zeta,\eta):=\sum_{C\in\zeta}\sum_{D\in\eta}\frac{|C\cap D|}{n}\cdot\log_{2}\left(\frac{|C \cap D| \cdot n}{|C|\cdot|D|}\right)$
	double MI = 0.0; // mutual information
//	const double nDouble = (double) n;
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
 					double factor1 =  ((double) sizeO) / totalOverlapDouble;
					assert ((size_old[C] * size_new[D]) != 0);
					DEBUG("overlap of " << C << " and " << D << " has size: " << sizeO);
					double frac2 = ((double) (sizeO * totalOverlap)) / (double) (size_old[C] * size_new[D]);
					assert (frac2 != 0);
					double factor2 = log_b(frac2, 2);
					DEBUG("frac2 = " << frac2 << ", factor1 = " << factor1 << ", factor2 = " << factor2);
					MI += factor1 * factor2;
				}
			}
		}
	}
#endif

	// sanity check
	assert (! std::isnan(MI));
	assert (MI >= 0.0);


	// compute entropy for both clusterings
	double H_old = entropy(oldClustering, totalOverlap, P_old);
	double H_new = entropy(newClustering, totalOverlap, P_new);

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
	DEBUG("sanity check, NMI: " << NMI);

	if (Aux::NumericTools::equal(NMID, 0.0)) {
		NMID = 0.0;
	}
	if (Aux::NumericTools::equal(NMID, 1.0)) {
		NMID = 1.0;
	}

	// if NMID is close to 0 because of numerical error
	if (! Aux::NumericTools::ge(NMID, 0.0)) {
		ERROR("Set NMID from below 0 to exactly 0: " << NMID);
		NMID = 0.0;
	}
	if (! Aux::NumericTools::le(NMID, 1.0)) {
		ERROR("Set NMID larger than 1 to exactly 1: " << NMID);
		NMID = 1.0;
	}

	assert (Aux::NumericTools::ge(NMID, 0.0));
	assert (Aux::NumericTools::le(NMID, 1.0));
}

std::vector<std::vector<count> > DynamicNMIDistance::confusionMatrix(Graph& G,
		Clustering& first, Clustering& second) {
	cluster firstUpperId = first.upperBound();
	cluster secondUpperId = second.upperBound();
	std::vector<std::vector<count> > confMatrix(firstUpperId);

	for (index i = 0; i < first.upperBound(); ++i) {
		confMatrix[i].resize(secondUpperId, 0);
	}

	TRACE("upperId in first, second: " << first.upperBound() << ", " << secondUpperId);

	second.forEntries([&](node u, cluster secondId) {
		if (this->isInBoth(u, first, second)) {
			TRACE("node " << u << ", id in first: " << first[u] << ", in second: " << second[u]);
			cluster firstId = first[u];
			assert(firstId < confMatrix.size() && secondId < confMatrix[firstId].size());
			confMatrix[firstId][secondId]++;
		}
	});

	return confMatrix;
}



} /* namespace NetworKit */
