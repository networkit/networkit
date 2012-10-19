/*
 * EdgeTripleGraphDataStructure.cpp
 *
 *  Created on: 16.10.2012
 *      Author: cls
 */

#include "log4cxx/logger.h"

#include "../globals.h"
#include "Graph.h"
#include "EdgeTripleGraphData.h"


namespace EnsembleClustering {

EdgeTuple::EdgeTuple() {
	this->i = 0;
	this->j = 0;
	this->w = 0.0;
}

EdgeTuple::EdgeTuple(id i, id j, double w) {
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

void EdgeTripleGraphData::connectNode(id v, std::vector<id> indices) {

	LOG4CXX_DEBUG(log4cxx::Logger::getRootLogger(), "Connecting node " << v << " with " << indices.size() << " other nodes");

	NodeTuple nodetuple;
	this->nodeData[v] = 0.0; // TODO: weight?

	int deg = indices.size();

	int i; // edge tuple array index
	EdgeTuple edgetuple;
	for (std::vector<id>::iterator iter = indices.begin(); iter != indices.end(); ++iter) {
		this->edgeData[i] = EdgeTuple(v, *iter, 0.0); // TODO: weight?
	}

}

} /* namespace EnsembleClustering */


