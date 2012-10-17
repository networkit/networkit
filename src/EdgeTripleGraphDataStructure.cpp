/*
 * EdgeTripleGraphDataStructure.cpp
 *
 *  Created on: 16.10.2012
 *      Author: cls
 */

#include "EdgeTripleGraphDataStructure.h"

namespace EnsembleClustering {

EdgeTuple::EdgeTuple(int i, int j, double w) {
	this->i = i;
	this->j = j;
	this->w = w;
}

NodeTuple::NodeTuple(double w) {
	this->w = w;
}

EdgeTripleGraphData::EdgeTripleGraphData() {
	// TODO Auto-generated constructor stub

}

EdgeTripleGraphData::~EdgeTripleGraphData() {
	// TODO Auto-generated destructor stub
}

} /* namespace EnsembleClustering */


