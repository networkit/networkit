/*
 * GraphData.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#include "GraphData.h"

namespace NetworKit {

namespace graph_impl {

count WeightedData::getMemoryUsage() const {
	count mem = 0;

	for (auto& w : edgeWeights) {
		mem += sizeof(edgeweight) * w.capacity();
	}

	return mem;
}

void WeightedData::shrinkToFit() {
	edgeWeights.shrink_to_fit();
	for (auto& w : edgeWeights) {
		w.shrink_to_fit();
	}
}

void WeightedData::addNode() {
	edgeWeights.push_back(std::vector<edgeweight>{});
}

count UndirectedData::getMemoryUsage() const {
	count mem = 0;
	mem += sizeof(count) * deg.capacity();
	for (auto& a : adja) {
		mem += sizeof(node) * a.capacity();
	}
	return mem;
}

void UndirectedData::shrinkToFit() {
	deg.shrink_to_fit();

	adja.shrink_to_fit();
	for (auto& a : adja) {
		a.shrink_to_fit();
	}
}

void UndirectedData::addNode() {
	deg.push_back(0);
	adja.push_back(std::vector<node>{});
}

index UndirectedData::indexInEdgeArray(node u, node v) const {
	for (index i = 0; i < adja[u].size(); i++) {
		node x = adja[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

count DirectedData::getMemoryUsage() const {
	count mem = 0;
	
	mem += sizeof(count) * inDeg.capacity();
	mem += sizeof(count) * outDeg.capacity();
	
	for (auto& a : inEdges) {
		mem += sizeof(node) * a.capacity();
	}
	for (auto& a : outEdges) {
		mem += sizeof(node) * a.capacity();
	}
	
	return mem;
}

void DirectedData::shrinkToFit() {
	inDeg.shrink_to_fit();
	outDeg.shrink_to_fit();
	
	inEdges.shrink_to_fit();
	for (auto& a : inEdges) {
		a.shrink_to_fit();
	}

	outEdges.shrink_to_fit();
	for (auto& a : outEdges) {
		a.shrink_to_fit();
	}
}

void DirectedData::addNode() {
	inDeg.push_back(0);
	outDeg.push_back(0);
	inEdges.push_back(std::vector<node>{});
	outEdges.push_back(std::vector<node>{});
}

index DirectedData::indexInInEdgeArray(node u, node v) const {
	for (index i = 0; i < inEdges[u].size(); i++) {
		node x = inEdges[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

index DirectedData::indexInOutEdgeArray(node u, node v) const {
	for (index i = 0; i < outEdges[u].size(); i++) {
		node x = outEdges[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

} /* namespace graph_impl */

} /* namespace NetworKit */