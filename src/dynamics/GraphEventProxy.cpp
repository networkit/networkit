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
	for (GraphEventHandler* handler : this->handlers) {
		handler->onNodeAddition(u);
	}
	return u;
}

void GraphEventProxy::removeNode(node u) {
	this->G->removeNode(u);
	for (GraphEventHandler* handler : this->handlers) {
		handler->onNodeRemoval(u);
	}
}

void GraphEventProxy::addEdge(node u, node v) {
	this->G->addEdge(u, v);
	for (GraphEventHandler* handler : this->handlers) {
		handler->onEdgeAddition(u, v);
	}
}

void GraphEventProxy::removeEdge(node u, node v) {
	this->G->removeEdge(u, v);
	for (GraphEventHandler* handler : this->handlers) {
		handler->onEdgeRemoval(u, v);
	}
}

void GraphEventProxy::setWeight(node u, node v, edgeweight w) {
	this->G->setWeight(u, v, w);
	for (GraphEventHandler* handler : this->handlers) {
		handler->onWeightUpdate(u, v, w);
	}
}

} /* namespace NetworKit */
