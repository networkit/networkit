/*
 * EdgeTripleGraphDataStructure.h
 *
 *  Created on: 16.10.2012
 *      Author: cls
 */

#ifndef EDGETRIPLEGRAPHDATASTRUCTURE_H_
#define EDGETRIPLEGRAPHDATASTRUCTURE_H_

#include "Graph.h"

namespace EnsembleClustering {




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
	double w; //!< edge weight

};


class NodeTuple {

public:

	NodeTuple();

	NodeTuple(double w);

	double w; //!< self-loop weight

};



class EdgeTripleGraphData {

public:

	EdgeTripleGraphData(int n, int m);

	virtual ~EdgeTripleGraphData();

	/**
	 * Connects node with to all its neighbors.
	 *
	 * @param[in]	v			node id
	 * @param[in]	indices		neighbor node ids
	 */
	virtual void connectNode(id v, std::vector<id> indices);

	int n; //!< number of nodes
	int m; //!< number of edges

	EdgeTuple* edgeData;
	NodeTuple* nodeData;


};

} /* namespace EnsembleClustering */
#endif /* EDGETRIPLEGRAPHDATASTRUCTURE_H_ */
