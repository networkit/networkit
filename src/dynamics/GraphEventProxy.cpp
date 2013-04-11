/*
 * GraphEventProxy.cpp
 *
 *  Created on: 03.04.2013
 *      Author: cls
 */

#include "GraphEventProxy.h"

namespace NetworKit {


GraphEventProxy::GraphEventProxy(Graph& G) {
	this->G = &G;
}

GraphEventProxy::~GraphEventProxy() {
	// TODO Auto-generated destructor stub
}

node GraphEventProxy::addNode() {
	node u = this->G->addNode();
	for (GraphEventHandler* observer : this->observers) {
		observer->onNodeAddition(u);
	}
	return u;
}

void GraphEventProxy::removeNode(node u) {
	this->G->removeNode(u);
	for (GraphEventHandler* observer : this->observers) {
		observer->onNodeRemoval(u);
	}
}

void GraphEventProxy::addEdge(node u, node v) {
	this->G->addEdge(u, v);
	for (GraphEventHandler* observer : this->observers) {
		observer->onEdgeAddition(u, v);
	}
}

void GraphEventProxy::removeEdge(node u, node v) {
	this->G->removeEdge(u, v);
	for (GraphEventHandler* observer : this->observers) {
		observer->onEdgeRemoval(u, v);
	}
}


void GraphEventProxy::setWeight(node u, node v, edgeweight w) {
	this->G->setWeight(u, v, w);
	for (GraphEventHandler* observer : this->observers) {
		observer->onWeightUpdate(u, v, w);
	}
}

void GraphEventProxy::registerObserver(GraphEventHandler* observer) {
	this->observers.push_back(observer);
}

} /* namespace NetworKit */
