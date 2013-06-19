/*
 * DynamicEnsemble.cpp
 *
 *  Created on: 19.06.2013
 *      Author: cls
 */

#include "DynamicEnsemble.h"

namespace NetworKit {

DynamicEnsemble::DynamicEnsemble() {
	// TODO Auto-generated constructor stub

}

DynamicEnsemble::~DynamicEnsemble() {
	// TODO Auto-generated destructor stub
}

void DynamicEnsemble::setGraph(Graph& G) {
}

Clustering DynamicEnsemble::run() {
}

std::string DynamicEnsemble::toString() const {
	return "DynamicEnsemble";
}

std::vector<count> DynamicEnsemble::getTimerHistory() {
}

void DynamicEnsemble::onNodeAddition(node u) {
}

void DynamicEnsemble::onNodeRemoval(node u) {
}

void DynamicEnsemble::onEdgeAddition(node u, node v, edgeweight w) {
}

void DynamicEnsemble::onEdgeRemoval(node u, node v, edgeweight w) {
}

void DynamicEnsemble::onWeightUpdate(node u, node v, edgeweight wOld,
		edgeweight wNew) {
}

void DynamicEnsemble::onTimeStep() {
}

} /* namespace NetworKit */
