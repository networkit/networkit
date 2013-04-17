/*
 * DynamicLabelPropagation.cpp
 *
 *  Created on: 27.03.2013
 *      Author: cls
 */

#include "DynamicLabelPropagation.h"

namespace NetworKit {

DynamicLabelPropagation::DynamicLabelPropagation(Graph& G, count theta) :
		DynamicClusterer(G),
		labels(G.numberOfNodes()),
		activeNodes(G.numberOfNodes()),
		weightedDegree(G.numberOfNodes(), 0.0),
		updateThreshold(theta),
		nUpdated(G.numberOfNodes()) {
	labels.allToSingletons(); // initialize labels to singleton clustering
	// PERFORMANCE: precompute and store incident edge weight for all nodes
	DEBUG("[BEGIN] Label Propagation: precomputing incident weight");
	this->G->parallelForNodes([&](node v) {
		weightedDegree[v] = this->G->weightedDegree(v);
	});
	DEBUG("[DONE] Label Propagation: precomputing incident weight");

}


DynamicLabelPropagation::~DynamicLabelPropagation() {
	// TODO Auto-generated destructor stub
}

void DynamicLabelPropagation::onNodeAddition(node u) {
	// update data structures
	activeNodes.push_back(true); // new node is active
	weightedDegree.push_back(0.0);
	labels.append(u);
	labels.toSingleton(u);
}

void DynamicLabelPropagation::onNodeRemoval(node u) {
	assert (G->degree(u) == 0);
	assert (weightedDegree[u] == 0.0);
	activeNodes[u] = false; //
}

void DynamicLabelPropagation::onEdgeAddition(node u, node v) {
	// update weighted degree
	edgeweight w = G->weight(u, v);
	weightedDegree[u] += w;
	weightedDegree[v] += w;
	// assert that this is consistent with the graph
	assert (G->weightedDegree(u) == weightedDegree[u]);
	assert (G->weightedDegree(v) == weightedDegree[v]);
	// activate source and target // TODO: strategy
	activeNodes[u] = true;
	activeNodes[v] = true;
}

void DynamicLabelPropagation::onEdgeRemoval(node u, node v) {
	// update weighted degree
	edgeweight w = G->weight(u, v);
	weightedDegree[u] -= w;
	weightedDegree[v] -= w;
	// assert that this is consistent with the graph
	assert (G->weightedDegree(u) == weightedDegree[u]);
	assert (G->weightedDegree(v) == weightedDegree[v]);
	// activate source and target // TODO: strategy
	this->activeNodes[u] = true;
	this->activeNodes[v] = true;
}


void DynamicLabelPropagation::onWeightUpdate(node u, node v, edgeweight wOld, edgeweight wNew) {
	//update weighted degree
	weightedDegree[u] += (wNew - wOld);
	weightedDegree[v] += (wNew - wOld);

	// activate source and target // TODO: strategy
	this->activeNodes[u] = true;
	this->activeNodes[v] = true;
}

std::string DynamicLabelPropagation::toString() const {
	return "DynamicLabelPropagation";
}

Clustering DynamicLabelPropagation::run() {

	Aux::Timer runtime;
	count nIterations = 0;
	runtime.start();
	while (nUpdated > updateThreshold) {
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
		// TODO:
	}

	runtime.stop();
	INFO("[DONE] LabelPropagation: iteration #" << nIterations << " - updated " << nUpdated << " labels, time spent: " << runtime.elapsedTag());

	return labels;
}



} /* namespace NetworKit */
