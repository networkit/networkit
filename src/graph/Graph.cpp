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
	TRACE("destructor of Graph " << this->getName() << " called: incorrect implementation leading to memory leak");
	// FIXME: stinger_free(this->stingerG);
}


Graph::Graph(const Graph& other) {
	this->n = other.n;
	this->name = other.name;
	this->stingerG = other.stingerG;	// FIXME: deep copy of STINGER data structure needed, which is difficult
	TRACE("copy: " << this->getName() << "(" << other.getName() << "): incorrect implementation possibly leading to shared state");

}

Graph& Graph::Graph::operator =(const Graph& other) {
	ERROR("assignment " << this->getName() << " = "<< other.getName() << ": operator= not yet implemented");
}



stinger* Graph::asSTINGER() const {
	return this->stingerG;
}

void Graph::insertEdge(node u, node v, double weight, int64_t type,
		int64_t timestamp) {
	// FIXME: do not store weight twice!
	// TODO: disable assertions for performance
	assert ((u <= this->n) && (v <= this->n));
	assert(u != v); // no self-loops allowed
	assert(!(this->hasEdge(u, v))); // no redundant edge insertions

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
	return stinger_total_edges(this->stingerG) / 2;
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
		#pragma omp atomic update
		total += this->weight(u, v);
	}, "parallel", "readonly");
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

double Graph::totalNodeWeight() {
	double total = 0.0;
	this->forallNodes([&](node v){
		total += this->weight(v);
	});
	return total;
}

double Graph::incidentWeight(node u) {
	double iw = 0.0;
	this->forallEdgesOf(u, [&](node u, node v){
		iw += this->weight(u, v);
	});
	return iw;
}

std::string Graph::toString() {
	std::stringstream strm;
	int64_t l = 0;	// number of weighted nodes (=self-loops)
	this->forallNodes([&](node v){
		if (this->weight(v) >= 0.0) {
			l += 1;
		}
	});
	strm << "Graph(name=" << this->getName() << ", n=" << this->numberOfNodes() << ", m=" << this->numberOfEdges() << ", l=" << l << ")";
	return strm.str();
}

} /* namespace EnsembleClustering */
