/*
 * Graph.h
 *
 *  Created on: 16.10.2012
 *      Author: cls
 */

#ifndef GRAPH_H_
#define GRAPH_H_

namespace EnsembleClustering {

typedef unsigned int id; //<! type definition of a node id


class Graph {


public:

	class Node;
	class Edge;

	Graph();

	virtual ~Graph();

	virtual int numberOfNodes();

	virtual int numberOfEdges();

	virtual int degree(Node v);

	virtual void neighbors(Node v);





};

} /* namespace EnsembleClustering */
#endif /* GRAPH_H_ */
