/*
 * EdgeTripeGraph.h
 *
 *  Created on: 31.10.2012
 *      Author: cls
 */

#ifndef EDGETRIPLEGRAPH_H_
#define EDGETRIPLEGRAPH_H_


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

	EdgeTuple(id i, id j, double w);

	id i; //!< edge source index
	id j; //!< edge target index
	weight w; //!< edge weight

};


class NodeTuple {

public:

	NodeTuple(weight w);

	weight w; //!< self-loop weight

};



class BucketTuple {



public:

	int begin;	//!< index into the EdgeTuple array pointing to the first EdgeTuple of this node
	int end;		//!< index into the EdgeTuple array pointing to the last EdgeTuple of this node

	BucketTuple(int begin, int end);


};

/** ------------------------------------------------------------------- **/

class EdgeTripleGraph: public EnsembleClustering::StaticGraphImplementor {

private:

	int n; 	//!< number of nodes
	int m; 	//!< number of edges

	std::vector<EdgeTuple> _edges;		//!< array containing EdgeTuples (ie. weighted edges grouped into buckets)
	std::vector<NodeTuple> _nodes;		//!< array containing NodeTuples (ie. self loop-weights)
	std::vector<BucketTuple> _buckets; 	//!< array containing BucketTuples (ie. begin and end indices into the edge array)

	// insert markers
	int ieh;	//!< insert edge here

public:

	EdgeTripleGraph(int n, int m);

	virtual ~EdgeTripleGraph();

	virtual int numberOfNodes();

	virtual int numberOfEdges();

	// virtual int degree(Node v) =0;

	// virtual void neighbors(Node v) =0;

	// virtual void nodes() =0;

	// virtual void edges() =0;

	// virtual void edges(Node u) =0;


	/****** non-interface methods *********/

	virtual void addAdjacencies(id u, std::vector<id> adjacencies);



};

} /* namespace EnsembleClustering */
#endif /* EDGETRIPLEGRAPH_H_ */
