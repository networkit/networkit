/*
 * Graph2.cpp
 *
 *  Created on: 04.02.2013
 *      Author: cls
 */

#include "Graph2.h"

namespace EnsembleClustering {

Graph2::Graph2() : defaultWeight(1.0) {
	// TODO: initialize
}

Graph2::~Graph2() {
	// TODO Auto-generated destructor stub
}


int64_t Graph2::find(node u, node v) const {
	int64_t i = none;
	for (node x : this->adja[u]) {
		i++;
		if (x == v) {
			return i;
		}
	}
	return none;
}

void Graph2::insertEdge(node u, node v) {
	this->adja[u].push_back(v);
	this->eweights[u].push_back(this->defaultWeight);
	// TODO: loop over all attributes, setting default attr
}

edgeweight Graph2::weight(node u, node v) const {
	int64_t vi = find(u, v);
	if (vi != none) {
		return this->eweights[u][vi];
	} else {
		return 0.0;
	}
}

void Graph2::setWeight(node u, node v, edgeweight w) {
	index vi = find(u, v);
	if (vi != none) {
		this->eweights[u][vi] = w;
	} else {
		// TODO: what if edge not there?
	}

}


bool Graph2::hasEdge(node u, node v) const {
	return (find(u, v) != -1);
}

node Graph2::addNode() {
	// TODO:
	// TODO: how to set capacity of std::vector
}

void Graph2::extendNodeRange(int64_t n) {
	// TODO:
}

bool Graph2::isEmpty() {
	return (n == 0);
}


int64_t Graph2::numberOfNodes() const {
	return this->n;
}


int64_t Graph2::numberOfEdges() const {
	// TODO:
}

} /* namespace EnsembleClustering */




