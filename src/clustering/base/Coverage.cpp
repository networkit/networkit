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



	double coverage = 0.0; // term $\frac{\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)}{\sum_{e \in E} \omega(e)}$
	double totalEdgeWeight = G.totalEdgeWeight(); // add edge weight

	if (totalEdgeWeight == 0.0) {
		throw std::invalid_argument(
				"Coverage is undefined for graphs without edges (including self-loops).");
	}

	IndexMap<cluster, double> intraEdgeWeight(zeta.upperBound(), 0.0); // cluster -> weight of its internal edges

	// compute intra-cluster edge weights per cluster
	// TODO: Make parallel, protect intraEdgeWeight[c]
	G.forEdgesAndTheirWeights(
			[&](node u, node v, edgeweight ew) {
				assert (u < zeta.numberOfNodes());
				assert (v < zeta.numberOfNodes());
				cluster c = zeta[u];
				cluster d = zeta[v];
				if (c == d) {
#ifdef DEBUG
					if ((c >= zeta.upperBound()) || (c < zeta.lowerBound())) {
						ERROR("c=" << c << " = zeta(" << u << ") is larger than upper bound: " << zeta.upperBound());
						ERROR("zeta: "); zeta.print();
					}
#endif
					assert ((zeta.lowerBound()) <= c && (c < zeta.upperBound()));
					intraEdgeWeight[c] += ew;
				} // else ignore edge
			});

	double intraEdgeWeightSum = 0.0; //!< term $\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)$
#pragma omp parallel for reduction(+:intraEdgeWeightSum)
	for (cluster c = zeta.lowerBound(); c < zeta.upperBound(); ++c) {
		intraEdgeWeightSum += intraEdgeWeight[c];
	}

	coverage = intraEdgeWeightSum / totalEdgeWeight;

	assert(coverage <= 1.0);
	assert(coverage >= 0.0);
	return coverage;
}

} /* namespace EnsembleClustering */
