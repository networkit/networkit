/*
 * Modularity.cpp
 *
 *  Created on: 10.12.2012
 *      Author: cls
 */

#include "Modularity.h"


namespace EnsembleClustering {


Modularity::Modularity(Graph& G) : QualityMeasure(G) {
	this->incidentWeight = new NodeMap<double>(this->G->numberOfNodes());
	this->precompute();
}

Modularity::~Modularity() {
	// TODO Auto-generated destructor stub
}

void Modularity::precompute() {
	int64_t n = G->numberOfNodes();
	// precompute incident edge weight for all nodes
	#pragma omp parallel for
	for (node v = 1; v <= n; ++v) {
		double iw = 0.0;
		READ_ONLY_FORALL_EDGES_OF_NODE_BEGIN((*G), v) {
			iw += G->weight(EDGE_SOURCE, EDGE_DEST);
		} READ_ONLY_FORALL_EDGES_OF_NODE_END();
		(*this->incidentWeight)[v] = iw;
	}
}

double Modularity::getQuality(const Clustering& zeta) {

	int64_t n = G->numberOfNodes();

	double coverage;	// term $\frac{\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)}{\sum_{e \in E} \omega(e)}$
	double expectedCoverage; // term $\frac{ \sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 }{4( \sum_{e \in E} \omega(e) )^2 }$
	double modularity; 	// mod = coverage - expectedCoverage

	double totalEdgeWeight = G->totalEdgeWeight();


	IndexMap<cluster, double> intraEdgeWeight(zeta.upperBound(), 0.0); // cluster -> weight of its internal edges


	G->forallEdges([&](node u, node v){
		cluster c = zeta[u];
		cluster d = zeta[v];
		if (c == d) {
			if (c > zeta.upperBound()) {
				DEBUG("c=" << c << ", upper bound=" << zeta.upperBound());
			}
			assert (c <= zeta.upperBound());
			#pragma omp atomic update
			intraEdgeWeight[c] += G->weight(u, v);
		} // else ignore edge
	}, "parallel", "readonly");


	double intraEdgeWeightSum = 0.0;	//!< term $\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)$
	for (cluster c = zeta.lowerBound(); c <= zeta.upperBound(); ++c) {
		intraEdgeWeightSum += intraEdgeWeight[c];
	}

	coverage = intraEdgeWeightSum / totalEdgeWeight;

	IndexMap<cluster, double> incidentWeightSum(zeta.upperBound(), 0.0);	//!< cluster -> sum of the weights of incident edges for all nodes


	for (node v = G->firstNode(); v <= n; ++v) {
		// add to cluster weight
		cluster c = zeta.clusterOf(v);
		incidentWeightSum[c] += (*this->incidentWeight)[v];
	}

	double totalIncidentWeight = 0.0; 	//!< term $\sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 $
	for (cluster c = zeta.lowerBound(); c <= zeta.upperBound(); ++c) {
		totalIncidentWeight += incidentWeightSum[c] * incidentWeightSum[c];	// squared
	}

	expectedCoverage = totalIncidentWeight / (4 * totalEdgeWeight * totalEdgeWeight);


	DEBUG("coverage: " << coverage);
	DEBUG("expected coverage: " << expectedCoverage);
	modularity = coverage - expectedCoverage;
	return modularity;
}

} /* namespace EnsembleClustering */
