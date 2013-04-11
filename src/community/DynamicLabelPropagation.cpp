/*
 * DynamicLabelPropagation.cpp
 *
 *  Created on: 27.03.2013
 *      Author: cls
 */

#include "DynamicLabelPropagation.h"

namespace NetworKit {

DynamicLabelPropagation::DynamicLabelPropagation(Graph& G) : DynamicClusterer(G) {
	this->G = &G;
	// this->labels(0);
}


DynamicLabelPropagation::~DynamicLabelPropagation() {
	// TODO Auto-generated destructor stub
}

void DynamicLabelPropagation::onNodeAddition(node u) {
	this->activeNodes.push_back(true); // new node is active
	// TODO: new node gets id as label
}

void DynamicLabelPropagation::onNodeRemoval(node u) {
	// TODO: implies removal of all incident edges
	this->G->forNeighborsOf(u, [&](node v){
		this->activeNodes[v] = true;
	});
	// TODO: node is no longer included in iteration
}

void DynamicLabelPropagation::onEdgeAddition(node u, node v) {
	this->activeNodes[u] = true;
	this->activeNodes[v] = true;
}

void DynamicLabelPropagation::onEdgeRemoval(node u, node v) {
	this->activeNodes[u] = true;
	this->activeNodes[v] = true;
}


void DynamicLabelPropagation::onWeightUpdate(node u, node v, edgeweight w) {
	this->activeNodes[u] = true;
	this->activeNodes[v] = true;
}

std::string DynamicLabelPropagation::toString() const {
	return "DynamicLabelPropagation";
}

Clustering DynamicLabelPropagation::run() {
	return Clustering(this->G->numberOfNodes()); // FIXME:

// TODO: standard PLP iteration

//	count nUpdated;
//
//	// propagate labels
//	while (nUpdated > this->updateThreshold) {
//		this->G->parallelForNodes([&](node u) {
//			// TODO:
//			});
//	}

}



} /* namespace NetworKit */
