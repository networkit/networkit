/*
 * BasicGraph.cpp
 *
 *  Created on: 20.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#include <stdexcept>
#include <sstream>

#include "BasicGraph.h"
#include "../auxiliary/Random.h"

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

template<Weighted w, Directed d>
BasicGraph<w, d>::BasicGraph(count n) :
	DData(n),
	n(n),
	m(0),
	z(n),
	t(0),
	exists(n, true) {
	
	// set name from global id
	static count nextGraphId = 1;
	std::stringstream sstm;
	sstm << "G#" << nextGraphId++;
	this->name = sstm.str();
}

template<>
std::string BasicGraph<Weighted::unweighted, Directed::undirected>::typ() const { return "Graph"; }

template<>
std::string BasicGraph<Weighted::weighted, Directed::undirected>::typ() const { return "WeightedGraph"; }

template<>
std::string BasicGraph<Weighted::unweighted, Directed::directed>::typ() const { return "DirectedGraph"; }

template<>
std::string BasicGraph<Weighted::weighted, Directed::directed>::typ() const { return "WeightedDirectedGraph"; }

template<Weighted w, Directed d>
count BasicGraph<w, d>::getMemoryUsage() const {

}

template<Weighted w, Directed d>
void BasicGraph<w, d>::shrinkToFit() {
	// ToDo implement
	throw std::runtime_error("TODO");
}


template<Weighted w, Directed d>
std::string BasicGraph<w, d>::toString() const {
	std::stringstream strm;
	strm << typ() << "(name=" << getName() << ", n=" << numberOfNodes()
			<< ", m=" << numberOfEdges() << ")";
	return strm.str();
}


/** NODE MODIFIERS **/

template<Weighted w, Directed d>
node BasicGraph<w, d>::addNode() {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
node BasicGraph<w, d>::addNode(float x, float y) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::removeNode(node v) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

/** NODE PROPERTIES **/

template<Weighted w, Directed d>
bool BasicGraph<w, d>::isIsolated(node v) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
edgeweight BasicGraph<w, d>::weightedDegree(node v) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
node BasicGraph<w, d>::randomNode() const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** EDGE MODIFIERS **/

template<Weighted w, Directed d>
void BasicGraph<w, d>::addEdge(node u, node v, edgeweight ew) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::removeEdge(node u, node v) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
bool BasicGraph<w, d>::hasEdge(node u, node v) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
node BasicGraph<w, d>::mergeEdge(node u, node v, bool discardSelfLoop) {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** GLOBAL PROPERTIES **/

template<Weighted w, Directed d>
count BasicGraph<w, d>::numberOfSelfLoops() const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** COORDINATES **/

template<Weighted w, Directed d>
void BasicGraph<w, d>::setCoordinate(node v, Point<float> value) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
Point<float>& BasicGraph<w, d>::getCoordinate(node v) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
float BasicGraph<w, d>::minCoordinate(count dim) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
float BasicGraph<w, d>::maxCoordinate(count dim) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::initCoordinates() {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** EDGE ATTRIBUTES **/

template<Weighted w, Directed d>
edgeweight BasicGraph<w, d>::weight(node u, node v) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::setWeight(node u, node v, edgeweight ew) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::increaseWeight(node u, node v, edgeweight ew) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
int BasicGraph<w, d>::addEdgeAttribute_double(double defaultValue) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
double BasicGraph<w, d>::attribute_double(node u, node v, int attrId) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::setAttribute_double(node u, node v, int attrId, double attr) {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** SUMS **/

template<Weighted w, Directed d>
edgeweight BasicGraph<w, d>::totalEdgeWeight() const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** Collections **/

template<Weighted w, Directed d>
std::vector<node> BasicGraph<w, d>::nodes() const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
std::vector<std::pair<node, node> > BasicGraph<w, d>::edges() const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


template<Weighted w, Directed d>
std::vector<node> BasicGraph<w, d>::neighbors(node u) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** NODE ITERATORS **/

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forNodes(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::parallelForNodes(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename C, typename L>
void BasicGraph<w, d>::forNodesWhile(C condition, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename C, typename L>
void BasicGraph<w, d>::forNodes(C condition, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forNodesInRandomOrder(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::balancedParallelForNodes(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forNodePairs(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::parallelForNodePairs(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** REDUCTION ITERATORS **/

template<Weighted w, Directed d>
template<typename L>
double BasicGraph<w, d>::parallelSumForNodes(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** EDGE ITERATORS **/

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forEdges(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::parallelForEdges(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forWeightedEdges(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::parallelForWeightedEdges(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forEdgesWithAttribute_double(int attrId, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** NEIGHBORHOOD ITERATORS **/

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forNeighborsOf(node u, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forWeightedNeighborsOf(node u, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forEdgesOf(node u, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forWeightedEdgesOf(node u, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** REDUCTION ITERATORS **/

template<Weighted w, Directed d>
template<typename L>
double BasicGraph<w, d>::parallelSumForWeightedEdges(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** GRAPH SEARCHES **/

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::BFSfrom(node r, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::BFSEdgesfrom(node r, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::DFSfrom(node r, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::DFSEdgesfrom(node r, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

} /* namespace graph_impl */

} /* namespace NetworKit */
