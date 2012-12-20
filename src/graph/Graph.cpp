/*
 * Graph.cpp
 *
 *  Created on: 28.11.2012
 *      Author: cls
 */

#include "Graph.h"

namespace EnsembleClustering {


Graph::Graph() {
	this->stingerG = stinger_new();
	this->nextNode = 1;
}

Graph::~Graph() {
}

Graph::Graph(stinger* stingerG) {
	this->stingerG = stingerG;
}



stinger* Graph::asSTINGER() const {
	return this->stingerG;
}

void Graph::insertEdge(node u, node v, double weight, int64_t type,
		int64_t timestamp) {
	// FIXME: do not store weight twice
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
	// TODO: is this sufficient? do isolated nodes have to be counted?
	// TODO: implement node counter
	return stinger_max_active_vertex(this->stingerG);
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
	node v = this->nextNode;
	this->nextNode += 1;
	// currently, Graph does not store of nodes any more than STINGER
	return v;
}

node Graph::lastNode() const {
	return this->numberOfNodes();
}

} /* namespace EnsembleClustering */
