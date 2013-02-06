/*
 * Graph.cpp
 *
 *  Created on: 04.02.2013
 *      Author: cls
 */

#include "Graph.h"

namespace EnsembleClustering {

Graph::Graph(count n) : n(n), maxn(n), deg(n, 0), adja(n), eweights(n) {
	// set name from global id
	static int64_t graphId = 1;
	std::stringstream sstm;
	sstm << "G#" << graphId++;
	this->name = sstm.str();
}

Graph::~Graph() {
	// TODO Auto-generated destructor stub
}


index Graph::find(node u, node v) const {
	index vi = none;
	for (node x : this->adja[u]) {
		vi++;
		if (x == v) {
			return vi;
		}
	}

	return none;
}

void Graph::insertEdge(node u, node v, edgeweight weight) {
	// set adjacency
	this->adja[u].push_back(v);
	this->adja[v].push_back(u);
	// increment degree counters
	this->deg[u] += 1;
	this->deg[v] += 1;
	// set edge weight
	this->eweights[u].push_back(weight);
	this->eweights[v].push_back(weight);
	// TODO: loop over all attributes, setving default attr
}

void Graph::removeEdge(node u, node v) {
	// remove adjacency
	index vi = find(u, v);
	index ui = find(v, u);
	if (vi == none) {
		ERROR("edge (" << u << "," << v << ") does not exist");
		// TODO: what if edge does not exist?
	} else {
		this->adja[u][vi] = none;
		this->adja[v][ui] = none; //FIXME:  assumpvion: u is at same index w.r.t. v as v w.r.t. u -
		// decrement degree counters
		this->deg[u] -= 1;
		this->deg[v] -= 1;
		// remove edge weight
		this->eweights[u][vi] = this->nullWeight;
		this->eweights[v][ui] = this->nullWeight;
		// TODO: remove attributes

	}
}

edgeweight Graph::weight(node u, node v) const {
	index vi = find(u, v);
	if (vi != none) {
		return this->eweights[u][vi];
	} else {
		return 0.0;
	}
}

void Graph::setWeight(node u, node v, edgeweight w) {
	index vi = find(u, v);
	if (vi != none) {
		this->eweights[u][vi] = w;
	} else {
		// TODO: what if edge not there?
	}
}

bool Graph::hasEdge(node u, node v) const {
	TRACE("find(" << u << "," << v << ") = " << find(u, v));
	return (find(u, v) != none);
}

node Graph::addNode() {
	this->n += 1;

	//update per node data structures
	this->deg.push_back(0);

	// update per edge data structures
	std::vector<node> adjacencyVector;	// vector of adjacencies for new node
	std::vector<edgeweight> edgeWeightVector;	// vector of edge weights for new node
	this->adja.push_back(adjacencyVector);
	this->eweights.push_back(edgeWeightVector);
}

void Graph::extendNodeRange(int64_t n) {
	throw std::runtime_error("TODO");
	// TODO:
}

bool Graph::isEmpty() {
	return (n == 0);
}

int64_t Graph::numberOfNodes() const {
	return this->n;
}

count Graph::degree(node v) const {
	return deg[v];
}

edgeweight Graph::weightedDegree(node v) const {
	throw std::runtime_error("DEPRECATED");
	// TODO:
}



int64_t Graph::numberOfEdges() const {
	// sum over all stored degrees
	// TODO: parallel sum?
	count mm = 0;
	this->forNodes([&](node v) {
		mm += this->deg[v];
	});
	count m = mm / 2;
	return m;
}

edgeweight Graph::totalEdgeWeight() {
	// TODO: optimize - replace by efficient parallel reduction?
	edgeweight sum = 0.0;
	this->forEdges([&](node u, node v){
		sum += this->weight(u, v);
	});
}

void Graph::setName(std::string name) {
	this->name = name;
}

std::string Graph::toString() {
	return "TODO";
	// TODO:
}

std::string Graph::getName() {
	// TODO: unneccessary if name becomes public attribute
	return this->name;
}

edgeweight Graph::totalNodeWeight() {
	throw std::runtime_error("DEPRECATED");
}

void Graph::setWeight(node u, edgeweight w) {
	WARN("DEPRECATED");
	this->setWeight(u, u, w);
}

edgeweight Graph::weight(node v) {
	WARN("DEPRECATED");
	return this->weight(v, v); // return self-loop weight
}

} /* namespace EnsembleClustering */




