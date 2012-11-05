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


public:

	class Node;
	class Edge;

	Graph();

	virtual ~Graph();

	virtual int numberOfNodes() =0;

	virtual int numberOfEdges() =0;

	virtual int degree(Node v) =0;

	virtual void neighbors(Node v) =0;





};

} /* namespace EnsembleClustering */
#endif /* GRAPH_H_ */
