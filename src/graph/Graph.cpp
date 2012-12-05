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
}

Graph::~Graph() {
}

Graph::Graph(stinger* stingerG) {
	this->stingerG = stingerG;
}



stinger* Graph::asSTINGER() {
	return this->stingerG;
}

void Graph::insertEdge(node u, node v, double weight, int64_t type,
		int64_t timestamp) {
	stinger_insert_edge_pair(this->stingerG, type, u, v, weight, timestamp);
}

bool Graph::hasEdge(node u, node v) {
	int to = stinger_has_typed_successor(this->stingerG, this->defaultEdgeType, u, v);
	int back = stinger_has_typed_successor(this->stingerG, this->defaultEdgeType, v, u);
	return to && back;
}


double Graph::getWeight(node v) {
	return stinger_vweight(this->stingerG, v);
}

double Graph::getWeight(edge uv) {
	return stinger_edgeweight(this->stingerG, uv.first, uv.second,
			this->defaultEdgeType);

}




double Graph::getWeight(node u, node v) {
	return stinger_edgeweight(this->stingerG, u, v, this->defaultEdgeType);
}

int64_t Graph::numberOfEdges() {
	return stinger_total_edges(this->stingerG);
}


int64_t Graph::numberOfNodes() {
	// TODO: is this sufficient? do isolated nodes have to be counted?
	return stinger_max_active_vertex(this->stingerG);
}

node Graph::firstNode() {
	return 1;
}

node Graph::lastNode() {
	return this->numberOfNodes();
}

} /* namespace EnsembleClustering */
