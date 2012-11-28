/*
 * Graph.cpp
 *
 *  Created on: 28.11.2012
 *      Author: cls
 */

#include "Graph.h"

namespace EnsembleClustering {


Graph::Graph(stinger* stingerG) {
	this->stingerG = stingerG;
}

Graph::~Graph() {
	// TODO Auto-generated destructor stub
}

stinger* Graph::asSTINGER() {
	return this->stingerG;
}

void Graph::insertEdge(node u, node v, double weight, int64_t type, int64_t timestamp) {
	stinger_insert_edge_pair(this->stingerG, type, u, v, weight, timestamp);
}

double Graph::getWeight(node v) {
	return stinger_vweight(this->stingerG, v);
}

double Graph::getWeight(edge uv) {
	return stinger_edgeweight(this->stingerG, uv.first, uv.second, this->defaultEdgeType);

}

double Graph::getWeight(node u, node v) {
	return stinger_edgeweight(this->stingerG, u, v, this->defaultEdgeType);
}

} /* namespace EnsembleClustering */
