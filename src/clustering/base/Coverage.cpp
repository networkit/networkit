/*
 * Coverage.cpp
 *
 *  Created on: 02.02.2013
 *      Author: cls
 */

#include "Coverage.h"

namespace EnsembleClustering {

Coverage::Coverage() {
	// TODO Auto-generated constructor stub

}

Coverage::~Coverage() {
	// TODO Auto-generated destructor stub
}

double Coverage::getQuality(const Clustering& zeta, Graph& G) {

	int64_t n = G.numberOfNodes();
	int64_t m = G.numberOfEdges();
	if ((m == 0) && (G.totalNodeWeight() == 0.0)) {
		throw std::invalid_argument("Coverage is undefined for graphs without edges (including self-loops).");
	}

	double coverage;	// term $\frac{\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)}{\sum_{e \in E} \omega(e)}$

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
				ERROR("zeta: "); zeta.print();
			}
			assert (c <= zeta.upperBound());
			intraEdgeWeight[c] += G.weight(u, v);
		} // else ignore edge
	});

	// .... and also add self-loops
	G.forallNodes([&](node v){
		cluster c = zeta.clusterOf(v);
		#pragma omp atomic update
		intraEdgeWeight[c] += G.weight(v);
	}, "parallel");


	double intraEdgeWeightSum = 0.0;	//!< term $\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)$
	for (cluster c = zeta.lowerBound(); c <= zeta.upperBound(); ++c) {
		intraEdgeWeightSum += intraEdgeWeight[c];
	}

	coverage = intraEdgeWeightSum / totalEdgeWeight;
	assert (coverage <= 1.0);
	assert (coverage >= 0.0);
	return coverage;
}

} /* namespace EnsembleClustering */
