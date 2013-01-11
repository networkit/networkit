/*
 * GraphFromAdjacencies.h
 *
 *  Created on: 27.11.2012
 *      Author: cls
 */

#ifndef GRAPHFROMADJACENCIES_H_
#define GRAPHFROMADJACENCIES_H_

#include <vector>
#include <utility>

#include "../aux/Log.h"
#include "../graph/Graph.h"

namespace EnsembleClustering {

/**
 * A 'builder' which constructs a STINGER-based graph from adjacencies.
 * An adjacency is a collection of node ids which represent a new node
 * as well as its incident edges.
 */
class GraphFromAdjacencies {

public:

	GraphFromAdjacencies();

	virtual ~GraphFromAdjacencies();

	/**
	 *
	 */
	virtual void readHeader(std::pair<int64_t, int64_t> header);

	/**
	 * Add next node and its adjacent edges.
	 */
	virtual void addAdjacencies(std::vector<node> adj);


	virtual Graph getGraph();

protected:

	Graph G;
	// int etype = 0;
	// double eweight = 1.0;
	// int timestamp = 0;
	node currentNode;

};


} /* namespace EnsembleClustering */
#endif /* STINGERFROMADJACENCIES_H_ */
