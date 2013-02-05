/*
 * Graph2.cpp
 *
 *  Created on: 04.02.2013
 *      Author: cls
 */

#include "Graph2.h"

namespace EnsembleClustering {

Graph2::Graph2(count n) : n(n), maxn(n), deg(n, 0), adja(n), eweights(n) {
}

Graph2::~Graph2() {
	// TODO Auto-generated destructor stub
}


index Graph2::find(node u, node v) const {
	index vi = none; // edge target index
	// NOTE: this is not symmetric
	for (node x : this->adja[u]) {
		vi++;
		if (x == v) {
			return vi;
		}
	}
	return none;
}

void Graph2::insertEdge(node u, node v) {
	// set adjacency
	this->adja[u].push_back(v);
	this->adja[v].push_back(u);
	// increment degree counters
	this->deg[u] += 1;
	this->deg[v] += 1;
	// set edge weight
	this->eweights[u].push_back(this->defaultEdgeWeight);
	this->eweights[v].push_back(this->defaultEdgeWeight);
	// TODO: loop over all attributes, setving default attr
}

void Graph2::removeEdge(node u, node v) {
	std::cout << std::endl;

	// remove adjacency
	index vi = find(u, v);
	index ui = find(v, u);
	if (vi == none) {
		ERROR("edge (" << u << "," << v << ") does not exist");
		// TODO: what if edge does not exist?
	} else {
		this->adja[u][vi] = none;
		this->adja[v][ui] = none;	//FIXME:  assumpvion: u is at same index w.r.t. v as v w.r.t. u -
		// decrement degree counters
		this->deg[u] -= 1;
		this->deg[v] -= 1;
		// remove edge weight
		this->eweights[u][vi] = this->nullWeight;
		this->eweights[v][ui] = this->nullWeight;
		// TODO: remove attributes

	}
}

edgeweight Graph2::weight(node u, node v) const {
	index vi = find(u, v);
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
	TRACE("find(" << u << "," << v << ") = " << find(u,v));
	return (find(u, v) != none);
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
	// sum over all stored degrees
	// TODO: parallel sum?
	count mm = 0;
	this->forNodes([&](node v) {
		mm += this->deg[v];
	});
	count m = mm / 2;
	return m;
}

} /* namespace EnsembleClustering */




