/*
 * ConcurrentGraph.cpp
 *
 *  Created on: 24.10.2013
 *      Author: cls
 */

#include "ConcurrentGraph.h"

namespace NetworKit {

ConcurrentGraph::ConcurrentGraph() : n(0), m(0), z(n), deg(z, 0), adja(z), eweights(z) {
}

ConcurrentGraph::ConcurrentGraph(count n) : n(n), m(0), z(n), deg(z, 0), adja(z), eweights(z) {
}

ConcurrentGraph::~ConcurrentGraph() {
	// TODO Auto-generated destructor stub
}

void ConcurrentGraph::addEdge(node u, node v, edgeweight weight) {
	assert (u >= 0);
	assert (u <= this->z);
	assert (v >= 0);
	assert (v <= this->z); // node ids must be in range
	if (u == v) { // self-loop case
		this->adja[u].push_back(u);
		this->deg[u] += 1;
		this->eweights[u].push_back(weight);
	} else {
		// set adjacency
		this->adja[u].push_back(v);
		this->adja[v].push_back(u);
		// increment degree counters
		this->deg[u] += 1;
		this->deg[v] += 1;
		// set edge weight
		this->eweights[u].push_back(weight);
		this->eweights[v].push_back(weight);
	}

	#pragma omp atomic update
	m++; // increasing the number of edges
}

void ConcurrentGraph::removeEdge(node u, node v) {
	// remove adjacency
	index ui = find(v, u);
	index vi = find(u, v);

	if (vi == none) {
		throw std::runtime_error("edge does not exist");
		// TODO: what if edge does not exist?
	} else {
		this->adja[u][vi] = none;
		this->adja[v][ui] = none;
		// decrement degree counters
		this->deg[u] -= 1;
		if (u != v) { // self-loops are counted only once
			this->deg[v] -= 1;
		}
		// remove edge weight
		this->eweights[u][vi] = this->nullWeight;
		this->eweights[v][ui] = this->nullWeight;

		#pragma omp atomic update
		m--; // decreasing the number of edges

	}
}

count ConcurrentGraph::numberOfEdges() const {
	return m;
}

count ConcurrentGraph::degree(node v) const {
	assert (v >= 0);
	assert (v <= z);
	return deg[v];
}

index ConcurrentGraph::find(node u, node v) const {
	for (index vi = 0; vi < this->adja[u].size(); ++vi) {
		node x = this->adja[u][vi];
		if (x == v) {
			return vi;
		}
	}
	return none;
}

} /* namespace NetworKit */
