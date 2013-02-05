/*
 * Modularity.cpp
 *
 *  Created on: 10.12.2012
 *      Author: cls
 */

#include "Modularity.h"


namespace EnsembleClustering {


Modularity::Modularity() : QualityMeasure() {
}

Modularity::~Modularity() {
	// TODO Auto-generated destructor stub
}


double Modularity::getQuality(const Clustering& zeta, Graph& G) {
	assert (G.numberOfNodes() <= zeta.numberOfNodes());

	int64_t n = G.numberOfNodes();
	int64_t m = G.numberOfEdges();
	if ((m == 0) && (G.totalNodeWeight() == 0.0)) {
		throw std::invalid_argument("Modularity is undefined for graphs without edges (including self-loops).");
	}

	double coverage;	// term $\frac{\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)}{\sum_{e \in E} \omega(e)}$
	double expectedCoverage; // term $\frac{ \sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 }{4( \sum_{e \in E} \omega(e) )^2 }$
	double modularity; 	// mod = coverage - expectedCoverage

	double totalEdgeWeight = 0.0;
	totalEdgeWeight += G.totalEdgeWeight(); // add edge weight
	totalEdgeWeight += G.totalNodeWeight();	// add "self-loop" weight


	NodeMap<double> incidentWeight(n, 0.0);
	// compute incident edge weight for all nodes
	G.forallNodes([&](node v){
		double iw = 0.0;
		G.forallEdgesOf(v, [&](node v, node w) {
			iw += G.weight(v, w);
		});
		iw += 2 * G.weight(v);	// Graph datastructure does not support self-loops. Node weights are used instead.
		// TODO: check DIMACS rules regarding self-loops
		incidentWeight[v] = iw;
	}, "parallel");


	IndexMap<cluster, double> intraEdgeWeight(zeta.upperBound(), 0.0); // cluster -> weight of its internal edges


	// compute intra-cluster edge weights per cluster
	G.forallEdges([&](node u, node v){
		assert (u <= zeta.numberOfNodes());
		assert (v <= zeta.numberOfNodes());
		cluster c = zeta[u];
		cluster d = zeta[v];
		if (c == d) {
			if (c > zeta.upperBound()) {
				ERROR("c=" << c << " = zeta(" << u << ") is larger than upper bound: " << zeta.upperBound());
				// ERROR("zeta: "); zeta.print();
			}
			assert (c <= zeta.upperBound());
			intraEdgeWeight[c] += G.weight(u, v);
		} // else ignore edge
	});

	// .... and also add self-loops
	G.forallNodes([&](node v){
		cluster c = zeta.clusterOf(v);
		intraEdgeWeight[c] += G.weight(v);
	});


	double intraEdgeWeightSum = 0.0;	//!< term $\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)$
	#pragma omp parallel for reduction(+:intraEdgeWeightSum)
	for (cluster c = zeta.lowerBound(); c <= zeta.upperBound(); ++c) {
		intraEdgeWeightSum += intraEdgeWeight[c];
	}

	coverage = intraEdgeWeightSum / totalEdgeWeight;

	IndexMap<cluster, double> incidentWeightSum(zeta.upperBound(), 0.0);	//!< cluster -> sum of the weights of incident edges for all nodes


	G.forallNodes([&](node v){
		// add to cluster weight
		cluster c = zeta[v];
		assert (zeta.lowerBound() <= c);
		assert (c <= zeta.upperBound());
		incidentWeightSum[c] += incidentWeight[v];
	});

	double totalIncidentWeight = 0.0; 	//!< term $\sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 $
	#pragma omp parallel for reduction(+:totalIncidentWeight)
	for (cluster c = zeta.lowerBound(); c <= zeta.upperBound(); ++c) {
		totalIncidentWeight += incidentWeightSum[c] * incidentWeightSum[c];	// squared
	}
	double divisor = 4 * totalEdgeWeight * totalEdgeWeight;
	assert (divisor != 0);	// do not divide by 0
	expectedCoverage = totalIncidentWeight / divisor;


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
