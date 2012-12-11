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
	// precompute incident edge weight for all nodes
	for (node v = 1; v <= G->numberOfNodes(); ++v) {
		double iw = 0.0;
		FORALL_EDGES_OF_NODE_BEGIN((*G), v) {
			iw += G->weight(EDGE_SOURCE, EDGE_DEST);
		} FORALL_EDGES_OF_NODE_END();
		(*this->incidentWeight)[v] = iw;
	}
}

double Modularity::getQuality(Clustering& zeta) {
	double coverage;	//!< term $\frac{\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)}{\sum_{e \in E} \omega(e)}$
	double expectedCoverage; //!< term $\frac{ \sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 }{4( \sum_{e \in E} \omega(e) )^2 }$
	double modularity; 	// mod = coverage - expectedCoverage

	double totalEdgeWeight = G->totalEdgeWeight();
	IndexMap<cluster, double> intraEdgeWeight(zeta.lastCluster()); //!< cluster -> weight of its internal edges

	// TODO: parallel
	FORALL_EDGES_BEGIN((*G)) {
		node u = EDGE_SOURCE;
		node v = EDGE_DEST;
		cluster c = zeta.clusterOf(u);
		cluster d = zeta.clusterOf(v);
		if (c == d) {
			// TODO: make critical section atomic
			intraEdgeWeight[c] += G->weight(u, v);
		} // else ignore edge
	} FORALL_EDGES_END();

	double intraEdgeWeightSum = 0.0;	//!< term $\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)$
	for (cluster c = zeta.firstCluster(); c <= zeta.lastCluster(); ++c) {
		intraEdgeWeightSum += intraEdgeWeight[c];
	}

	coverage = intraEdgeWeightSum / totalEdgeWeight;

	IndexMap<cluster, double> incidentWeightSum(zeta.lastCluster());	//!< cluster -> sum of the weights of incident edges for all nodes


	for (node v = G->firstNode(); v <= G->numberOfNodes(); ++v) {
		// add to cluster weight
		cluster c = zeta.clusterOf(v);
		incidentWeightSum[c] += (*this->incidentWeight)[v];
	}

	double totalIncidentWeight = 0.0; 	//!< term $\sum_{C \in \zeta}( \sum_{v \in C} \omega(v) )^2 $
	for (cluster c = zeta.firstCluster(); c <= zeta.lastCluster(); ++c) {
		totalIncidentWeight += incidentWeightSum[c] * incidentWeightSum[c];	// squared
	}

	expectedCoverage = totalIncidentWeight / (4 * totalEdgeWeight * totalEdgeWeight);

	modularity = coverage - expectedCoverage;
	return modularity;
}

} /* namespace EnsembleClustering */
