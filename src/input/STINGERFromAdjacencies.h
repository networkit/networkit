/*
 * STINGERFromAdjacencies.h
 *
 *  Created on: 27.11.2012
 *      Author: cls
 */

#ifndef STINGERFROMADJACENCIES_H_
#define STINGERFROMADJACENCIES_H_

#include <vector>

extern "C" {
	#include "stinger.h"
}

typedef stinger graph;
typedef int64_t node;

namespace EnsembleClustering {

class STINGERFromAdjacencies {

public:

	STINGERFromAdjacencies();

	virtual ~STINGERFromAdjacencies();

	/**
	 * Create new STINGER instance.
	 */
	virtual void createGraph();

	/**
	 * Add next node and its adjacent edges.
	 */
	virtual void addAdjacencies(std::vector<node> adj);

	virtual graph* getGraph();

protected:

	graph* G;
	int etype = 0;
	double eweight = 1.0;
	int timestamp = 0;
	node currentNode;

};


} /* namespace EnsembleClustering */
#endif /* STINGERFROMADJACENCIES_H_ */
