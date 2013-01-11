/*
 * Graph.cpp
 *
 *  Created on: 28.11.2012
 *      Author: cls
 */

#include "Graph.h"

namespace EnsembleClustering {


Graph::Graph(int64_t n) {
	this->stingerG = stinger_new(); // TODO: manage memory
	this->n = n;
	static int64_t graphId = 1;
	std::stringstream sstm;
	sstm << "noname#" << graphId++;
	this->name = sstm.str();
}

Graph::~Graph() {
	// TODO: destructor stub
}



stinger* Graph::asSTINGER() const {
	return this->stingerG;
}

void Graph::insertEdge(node u, node v, double weight, int64_t type,
		int64_t timestamp) {
	// FIXME: do not store weight twice (?)
	assert ((u <= this->n) && (v <= this->n));	// TODO: disable assertions for performance
	assert(u != v); // no self-loops allowed
	stinger_insert_edge_pair(this->stingerG, type, u, v, weight, timestamp);
}

bool Graph::hasEdge(node u, node v) const {
	int to = stinger_has_typed_successor(this->stingerG, this->defaultEdgeType, u, v);
	int back = stinger_has_typed_successor(this->stingerG, this->defaultEdgeType, v, u);
	return to && back;
}


double Graph::weight(node v) const {
	return stinger_vweight(this->stingerG, v);
}

double Graph::weight(edge uv) const {
	return stinger_edgeweight(this->stingerG, uv.first, uv.second,
			this->defaultEdgeType);

}


int64_t Graph::numberOfEdges() const {
	return stinger_total_edges(this->stingerG);
}


int64_t Graph::numberOfNodes() const {
	return this->n;
}

node Graph::firstNode() const {
	return 1;
}

int64_t Graph::degree(node u) const {
	int64_t deg = stinger_outdegree(this->stingerG, u);
	// each ndirected edge is represented by two directed edges
	assert (deg == stinger_indegree(this->stingerG, u));
	return deg;
}

double Graph::totalEdgeWeight() {
	double total = 0.0;
	this->forallEdges([&](node u, node v) {
		total += this->weight(u, v);
	}, "readonly");
	return total;
}

void Graph::removeEdge(node u, node v) {
	stinger_remove_edge_pair(this->stingerG, this->defaultEdgeType, u, v);
}

node Graph::addNode() {
	this->n += 1;	// increment node range by one
	node v = this->n;
	return v;
}

void Graph::extendNodeRange(int64_t n) {
	assert(n > this->n); // new range must be greater than old range
	this->n = n;
}

bool Graph::isEmpty() {
	bool empty = (this->n == 0);
	assert(this->numberOfEdges() == 0);	// stinger edge data structure should not contain any edges either
	return empty;
}

node Graph::lastNode() const {
	return this->numberOfNodes();
}

void Graph::setName(std::string name) {
	this->name = name;
}

std::string Graph::getName() const {
	return this->name;
}

std::string Graph::toString() const {
	std::stringstream strm;
	strm << "Graph(name=" << this->getName() << ", n=" << this->numberOfNodes() << ", m=" << this->numberOfEdges() << ")";
	return strm.str();
}

} /* namespace EnsembleClustering */
