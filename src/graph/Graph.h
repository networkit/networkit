/*
 * Graph.h
 *
 *  Created on: 16.10.2012
 *      Author: cls
 */

#ifndef GRAPH_H_
#define GRAPH_H_

namespace EnsembleClustering {



class Graph {

protected:

	GraphImplementor* implementor;


public:

	class Node;
	class Edge;

	// TODO: templating instead of typedef weight?

	Graph();

	virtual ~Graph();

	virtual int numberOfNodes() =0;

	virtual int numberOfEdges() =0;

	virtual int degree(Node v) =0;

	virtual void neighbors(Node v) =0;

	virtual void nodes() =0;

	virtual void edges() =0;

	virtual void edges(Node u) =0;


};



} /* namespace EnsembleClustering */
#endif /* GRAPH_H_ */
