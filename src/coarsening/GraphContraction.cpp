/*
 * GraphContraction.cpp
 *
 *  Created on: 25.12.2012
 *      Author: cls
 */

#include "GraphContraction.h"

namespace EnsembleClustering {

GraphContraction::GraphContraction(Graph& fine, Graph& coarse, NodeMap<node>& fineToCoarse) {
	this->fine = &fine;
	this->coarse = &coarse;
	this->fineToCoarse = &fineToCoarse;
}

GraphContraction::~GraphContraction() {
	// TODO Auto-generated destructor stub
}

Graph& GraphContraction::getFineGraph() {
	return *(this->fine);
}

Graph& GraphContraction::getCoarseGraph() {
	return *(this->coarse);
}

NodeMap<node>& GraphContraction::getFineToCoarseMap() {
	return *(this->fineToCoarse);
}

} /* namespace EnsembleClustering */
