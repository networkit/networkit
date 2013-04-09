/*
 * GraphEventProxy.h
 *
 *  Created on: 03.04.2013
 *      Author: cls
 */

#ifndef GRAPHEVENTPROXY_H_
#define GRAPHEVENTPROXY_H_

#include "../graph/Graph.h"
#include "GraphEventHandler.h"

namespace NetworKit {

class GraphEventProxy {

protected:

	std::vector<GraphEventHandler*> handlers;


public:

	Graph* G;

	GraphEventProxy(Graph& G);

	virtual ~GraphEventProxy();

	virtual node addNode();

	virtual void removeNode(node u);

	virtual void addEdge(node u, node v);

	virtual void removeEdge(node u, node v);

	virtual void setWeight(node u, node v, edgeweight w);
};

} /* namespace NetworKit */
#endif /* GRAPHEVENTPROXY_H_ */
