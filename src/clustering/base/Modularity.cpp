/*
 * Modularity.cpp
 *
 *  Created on: 10.12.2012
 *      Author: cls
 */

#include "Modularity.h"
#include "Coverage.h"


namespace EnsembleClustering {


Modularity::Modularity() : QualityMeasure() {
}

Modularity::~Modularity() {
	// TODO Auto-generated destructor stub
}


double Modularity::getQuality(const Clustering& zeta, Graph& G) {
	assert (G.numberOfNodes() <= zeta.numberOfNodes());



	Coverage cov;
	double coverage = cov.getQuality(zeta, G); // deprecated: intraEdgeWeightSum / totalEdgeWeight;
	double expectedCoverage; // term $\frac{ \sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 }{4( \sum_{e \in E} \omega(e) )^2 }$
	double modularity; 	// mod = coverage - expectedCoverage
	double totalEdgeWeight = G.totalEdgeWeight(); // add edge weight

	if (totalEdgeWeight == 0.0) {
		ERROR("G: m=" << G.numberOfEdges() << "n=" << G.numberOfNodes());
		throw std::invalid_argument("Modularity is undefined for graphs without edges (including self-loops).");
	}

	IndexMap<cluster, double> incidentWeightSum(zeta.upperBound(), 0.0);	//!< cluster -> sum of the weights of incident edges for all nodes

	// compute volume of each cluster
	G.forNodes([&](node v){
		// add to cluster weight
		cluster c = zeta[v];
		assert (zeta.lowerBound() <= c);
		assert (c < zeta.upperBound());
		incidentWeightSum[c] += G.weightedDegree(v) + G.weight(v,v); // account for self-loops a second time
	});

	// compute sum of squared cluster volumes and divide by squared graph volume
	// double totalIncidentWeight = 0.0; 	//!< term $\sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 $
	expectedCoverage = 0.0;
	double divisor = 4 * totalEdgeWeight * totalEdgeWeight;
	assert (divisor != 0);	// do not divide by 0

	#pragma omp parallel for reduction(+:expectedCoverage)
	for (cluster c = zeta.lowerBound(); c < zeta.upperBound(); ++c) {
		expectedCoverage += incidentWeightSum[c] * incidentWeightSum[c] / divisor;	// squared
	}

	TRACE("coverage: " << coverage);
	TRACE("expected coverage: " << expectedCoverage);

	// assert ranges of coverage
	assert(coverage <= 1.0);
	assert(coverage >= 0.0);
	assert(expectedCoverage <= 1.0);
	assert(expectedCoverage >= 0.0);

	modularity = coverage - expectedCoverage;

	assert(! std::isnan(modularity));	// do not return NaN
	// do not return anything not in the range of modularity values
	assert(modularity >= -0.5);
	assert(modularity <= 1);
	return modularity;
}

} /* namespace EnsembleClustering */
