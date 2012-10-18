/*
 * EdgeTripleGraphDataStructure.cpp
 *
 *  Created on: 16.10.2012
 *      Author: cls
 */

#include "EdgeTripleGraphData.h"


namespace EnsembleClustering {

EdgeTuple::EdgeTuple() {
	this->i = 0;
	this->j = 0;
	this->w = 0.0;
}

EdgeTuple::EdgeTuple(int i, int j, double w) {
	this->i = i;
	this->j = j;
	this->w = w;
}

NodeTuple::NodeTuple() {
	this->w = 0.0;
}

NodeTuple::NodeTuple(double w) {
	this->w = w;
}



/*** Graph Data Structure here ***/

EdgeTripleGraphData::EdgeTripleGraphData(int n, int m) {

	this->n = n;
	this->m = m;

	this->nodeData = new NodeTuple[n];
	this->edgeData = new EdgeTuple[m];
}

EdgeTripleGraphData::~EdgeTripleGraphData() {
	// TODO Auto-generated destructor stub
}

} /* namespace EnsembleClustering */


