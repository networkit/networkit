/*
 * EdgeTripeGraph.h
 *
 *  Created on: 31.10.2012
 *      Author: cls
 */

#ifndef EDGETRIPLEGRAPH_H_
#define EDGETRIPLEGRAPH_H_

#include "Graph.h"

#include <vector>

namespace EnsembleClustering {

typedef unsigned int id; //<! type definition of a node id
typedef double weight; //<! double weights

/**
 * EdgeTuple
 * 		contains triple (i, j, w) where
 * 			i: edge source index
 * 			j: edge target index
 * 			w: edge weight
 */
class EdgeTuple {

public:

	EdgeTuple();

	EdgeTuple(id i, id j, double w);

	id i; //!< edge source index
	id j; //!< edge target index
	weight w; //!< edge weight

};


class NodeTuple {

public:

	NodeTuple();

	NodeTuple(weight w);

	weight w; //!< self-loop weight

};



class BucketTuple {


	id begin;	//!< index into the EdgeTuple array pointing to the first EdgeTuple of this node
	id end;		//!< index into the EdgeTuple array pointing to the last EdgeTuple of this node

};


class EdgeTripleGraph: public EnsembleClustering::Graph {

private:

	std::vector<EdgeTuple> edges;		//!< array containing EdgeTuples (ie. weighted edges grouped into buckets)
	std::vector<NodeTuple> nodes;		//!< array containing NodeTuples (ie. self loop-weights)
	std::vector<BucketTuple> buckets; 	//!< array containing BucketTuples (ie. begin and end indices into the edge array)

public:
	EdgeTripleGraph();
	virtual ~EdgeTripleGraph();
};

} /* namespace EnsembleClustering */
#endif /* EDGETRIPLEGRAPH_H_ */
