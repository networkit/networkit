/*
 * EdgeTripleGraphDataStructure.h
 *
 *  Created on: 16.10.2012
 *      Author: cls
 */

#ifndef EDGETRIPLEGRAPHDATASTRUCTURE_H_
#define EDGETRIPLEGRAPHDATASTRUCTURE_H_

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

	EdgeTuple(int i, int j, double w);

	int i; //!< edge source index
	int j; //!< edge target index
	double w; //!< edge weight

};


class NodeTuple {

public:

	NodeTuple(double w);

	double w; //!< self-loop weight

};



class EdgeTripleGraphData {

public:

	EdgeTripleGraphData();

	virtual ~EdgeTripleGraphData();

	EdgeTuple edgeArray[];

	NodeTuple nodeArray[];


};

} /* namespace EnsembleClustering */
#endif /* EDGETRIPLEGRAPHDATASTRUCTURE_H_ */
