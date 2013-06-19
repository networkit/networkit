/*
 * DynamicLabelPropagation.cpp
 *
 *  Created on: 27.03.2013
 *      Author: cls
 */

#include "DynamicLabelPropagation.h"

namespace NetworKit {

DynamicLabelPropagation::DynamicLabelPropagation() {
#ifndef CORE
	throw std::runtime_error("Nullary constructor needed only for Python Shell - no proper initialization");
#endif
}

DynamicLabelPropagation::DynamicLabelPropagation(count theta, std::string strategyName) :
		DynamicCommunityDetector(),
		updateThreshold(theta),
		nUpdated(0) {

	this->G = NULL; // G is set in method setGraph
	// select prep strategy
	if (strategyName == "Reactivate") {
		this->prepStrategy = new DynamicLabelPropagation::Reactivate(this);
	} else if (strategyName == "ReactivateNeighbors") {
		this->prepStrategy = new DynamicLabelPropagation::ReactivateNeighbors(this);
	} else {
		throw std::runtime_error("unknown prep strategy");
	}

}


DynamicLabelPropagation::~DynamicLabelPropagation() {
	// TODO Auto-generated destructor stub
}

void DynamicLabelPropagation::onNodeAddition(node u) {
	// update data structures
	activeNodes.push_back(true); // new node is active
	weightedDegree.push_back(0.0);
	labels.append(u); // extend label array by 1 entry

	prepStrategy->onNodeAddition(u);
}

void DynamicLabelPropagation::onNodeRemoval(node u) {
	assert (G->degree(u) == 0);
	assert (weightedDegree[u] == 0.0);

	prepStrategy->onNodeRemoval(u);
}

void DynamicLabelPropagation::onEdgeAddition(node u, node v, edgeweight w) {
	// update weighted degree
	if (u != v) {
		weightedDegree[u] += w;
		weightedDegree[v] += w;
	} else {
		weightedDegree[u] += w; // self-loop case
	}

	// assert that this is consistent with the graph
	assert (G->weightedDegree(u) == weightedDegree[u]);
	assert (G->weightedDegree(v) == weightedDegree[v]);
	// apply prep strategy
	this->prepStrategy->onEdgeAddition(u, v);
}

void DynamicLabelPropagation::onEdgeRemoval(node u, node v, edgeweight w) {
	// update weighted degree
	if (u != v) {
		weightedDegree[u] -= w;
		weightedDegree[v] -= w;
	} else {
		weightedDegree[u] -= w; // self-loop case
	}
	// assert that this is consistent with the graph
	assert (G->weightedDegree(u) == weightedDegree[u]);
	assert (G->weightedDegree(v) == weightedDegree[v]);
	// apply prep strategy
	this->prepStrategy->onEdgeRemoval(u, v);
}


void DynamicLabelPropagation::onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew) {
	//update weighted degree
	if (u != v) {
		weightedDegree[u] += (wNew - wOld);
		weightedDegree[v] += (wNew - wOld);
	} else {
		weightedDegree[u] += (wNew - wOld);
	}

	this->prepStrategy->onWeightUpdate(u, v, wOld, wNew);
}

void DynamicLabelPropagation::setGraph(Graph& G) {
	if (!G.isEmpty()) {
		throw std::runtime_error("G is not an empty graph. Currently, it is assumed that this algorithm is initialized with an empty graph, which is then constructed incrementally");
	}

	this->G = &G;
}

void DynamicLabelPropagation::onTimeStep() {
	// ignore
}


std::string DynamicLabelPropagation::toString() const {
	std::stringstream strm;
	strm << "DynamicLabelPropagation(" << this->prepStrategy->toString() << ")";
	return strm.str();
}

Clustering DynamicLabelPropagation::run() {
	if (this->G == NULL) {
		throw std::runtime_error("pointer to current graph was not initialized - call setGraph first");
	}

	INFO("running DynamicLabelPropagation at t=" << G->time());

	Aux::Timer runtime;
	count nIterations = 0;
	nUpdated = G->numberOfNodes(); // starts while loop - TODO: correct?

	runtime.start();
	while (nUpdated > updateThreshold) {
		nUpdated = 0;
		nIterations += 1;
		G->parallelForNodes([&](node u){
			if ((activeNodes[u]) && (G->degree(u) > 0)) {
				std::map<label, double> labelWeights; // neighborLabelCounts maps label -> frequency in the neighbors

				// weigh the labels in the neighborhood of v
				G->forWeightedNeighborsOf(u, [&](node v, edgeweight weight) {
					label lv = labels[v];
					labelWeights[lv] += weight; // add weight of edge {v, w}
				});

				// get heaviest label
				label heaviest = std::max_element(labelWeights.begin(),
								labelWeights.end(),
								[](const std::pair<label, edgeweight>& p1, const std::pair<label, edgeweight>& p2) {
									return p1.second < p2.second;})->first;

				if (labels[u] != heaviest) { // UPDATE
					labels[u] = heaviest;
					nUpdated += 1; // TODO: atomic update?
					G->forNeighborsOf(u, [&](node v) {
						activeNodes[v] = true;
					});
				} else {
					activeNodes[u] = false;
				}
			} else { /* node is isolated */ }
		}); // end parallel for nodes
		TRACE("nUpdated = " << nUpdated);
	} // end while

	runtime.stop();
	INFO("[DONE] LabelPropagation: iteration #" << nIterations << " - updated " << nUpdated << " labels, time spent: " << runtime.elapsedTag());

	return labels;
}


// PREP STRATEGY IMPLEMENTATIONS

DynamicLabelPropagation::Reactivate::Reactivate(DynamicLabelPropagation* dynPLP) {
	this->dynPLP = dynPLP;
}

DynamicLabelPropagation::Reactivate::~Reactivate() {
}


void DynamicLabelPropagation::Reactivate::onNodeAddition(node u) {
	dynPLP->labels.toSingleton(u);
	dynPLP->activeNodes[u] = true;
	assert (dynPLP->G->degree(u) == 0); // new node has no incident edges
}

void DynamicLabelPropagation::Reactivate::onNodeRemoval(node u) {
	dynPLP->activeNodes[u] = false; // assumption: this node can never be reactivated
	// incident edges must have been removed before, so neighborhood is already active at this point
}

void DynamicLabelPropagation::Reactivate::onEdgeAddition(node u, node v, edgeweight w) {
	// activate the affected nodes
	dynPLP->activeNodes[u] = true;
	dynPLP->activeNodes[v] = true;
}

void DynamicLabelPropagation::Reactivate::onEdgeRemoval(node u, node v, edgeweight w) {
	this->onEdgeAddition(u, v, w);
}

void DynamicLabelPropagation::Reactivate::onWeightUpdate(node u, node v,
		edgeweight wOld, edgeweight wNew) {
	this->onEdgeAddition(u, v);
}


DynamicLabelPropagation::ReactivateNeighbors::ReactivateNeighbors(DynamicLabelPropagation* dynPLP) {
	this->dynPLP = dynPLP;
}

DynamicLabelPropagation::ReactivateNeighbors::~ReactivateNeighbors() {
}

void DynamicLabelPropagation::ReactivateNeighbors::onNodeAddition(node u) {
	dynPLP->labels.toSingleton(u);
	dynPLP->activeNodes[u] = true;
	assert (dynPLP->G->degree(u) == 0); // new node has no incident edges
}

void DynamicLabelPropagation::ReactivateNeighbors::onNodeRemoval(node u) {
	dynPLP->activeNodes[u] = false; // assumption: this node can never be reactivated
	// incident edges must have been removed before, so neighborhood is already active at this point
}

void DynamicLabelPropagation::ReactivateNeighbors::onEdgeAddition(node u, node v, edgeweight w) {
	// activate the affected nodes and their neighbors
	dynPLP->activeNodes[u] = true;
	dynPLP->activeNodes[v] = true;
	dynPLP->G->forNeighborsOf(u, [&](node x) {
		dynPLP->activeNodes[x] = true;
	});
	dynPLP->G->forNeighborsOf(v, [&](node x) {
		dynPLP->activeNodes[x] = true;
	});
}

void DynamicLabelPropagation::ReactivateNeighbors::onEdgeRemoval(node u, node v, edgeweight w) {
	this->onEdgeAddition(u, v, w);
}

void DynamicLabelPropagation::ReactivateNeighbors::onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew) {
	this->onEdgeAddition(u, v, wNew);
}

std::string DynamicLabelPropagation::Reactivate::toString() {
	return "Reactivate";
}

void DynamicLabelPropagation::Reactivate::onTimeStep() {
}

std::string DynamicLabelPropagation::ReactivateNeighbors::toString() {
	return "ReactivateNeighbors";
}

void DynamicLabelPropagation::ReactivateNeighbors::onTimeStep() {
}

} /* namespace NetworKit */


