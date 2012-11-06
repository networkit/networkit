/*
 * EdgeTripleGraph.cpp
 *
 *  Created on: 31.10.2012
 *      Author: cls
 */

#include "EdgeTripleGraph.h"

namespace EnsembleClustering {


EdgeTripleGraph::EdgeTripleGraph(int n, int m) {
	this->n = n;
	this->m = m;
	// make room for data tuples
	this->_nodes.resize(n);
	this->_buckets.resize(n);
	this->_edges.resize(m);
	// set insert markers
	this->ieh = 0;

}

EdgeTripleGraph::~EdgeTripleGraph() {
}

int EdgeTripleGraph::numberOfNodes() {
	return this->n;
}



int EdgeTripleGraph::numberOfEdges() {
	return this->m;
}

void EdgeTripleGraph::addAdjacencies(id u, std::vector<id> adjacencies) {

	// buckets
	int deg = adjacencies.size();
	this->_buckets[u] = new EnsembleClustering::BucketTuple(this->ieh, this->ieh + deg);
	this->ieh += deg;

	// nodes
	double w = 0.0; // FIXME: allow weighted self-loops
	this->_nodes[u] = new EnsembleClustering::NodeTuple(w);

	// edges
	id v;
	for (int i = 0; i < deg; ++i) {
		v = adjacencies[i];
		// if ids both even or odd
		if (u % 2 == v % 2) {
			// TODO:
		} else {

		}
	}

}

} /* namespace EnsembleClustering */

EnsembleClustering::BucketTuple::BucketTuple(int begin, int end) {
}
