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


/** PRIVATE HELPERS **/

template<Weighted w>
index find_impl(const BasicGraph<w, Directed::undirected>& G, node u, node v) {
	for (index i = 0; i < G.adja[u].size(); i++) {
		node x = G.adja[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

template<Weighted w>
index find_impl(const BasicGraph<w, Directed::directed>& G, node u, node v) {
	for (index i = 0; i < G.outEdges[u].size(); i++) {
		node x = G.outEdges[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}


/** GRAPH INFORMATION **/

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
	count mem = WData::getMemoryUsage() + DData::getMemoryUsage();

	mem += exists.capacity() / 8;

	for (auto& map : edgeMaps_double) {
		for (auto& a : map) {
			mem += sizeof(double) * a.capacity();
		}
	}

	return mem;
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::shrinkToFit() {
	WData::shrinkToFit();
	DData::shrinkToFit();

	exists.shrink_to_fit();

	edgeMaps_double.shrink_to_fit();
	for (auto& map : edgeMaps_double) {
		map.shrink_to_fit();
		for (auto& a : map) {
			a.shrink_to_fit();
		}
	}

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

template<Weighted w>
bool isIsolated_impl(const BasicGraph<w, Directed::undirected>& G, node v) {
	return G.deg[v] == 0;
}

template<Weighted w>
bool isIsolated_impl(const BasicGraph<w, Directed::directed>& G, node v) {
	return G.inDeg[v] == 0 && G.outDeg[v];
}

template<Weighted w, Directed d>
edgeweight BasicGraph<w, d>::weightedDegree(node v) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
node BasicGraph<w, d>::randomNode() const {
	assert (numberOfNodes() > 0);

	node v;
	do {
		v = Aux::Random::integer(z);
	} while (!exists[v]);

	return v;
}


/** EDGE MODIFIERS **/

template<Weighted w>
void addEdge_impl(BasicGraph<w, Directed::undirected>& G, node u, node v, edgeweight ew) {
	assert (u >= 0);
	assert (u < G.z);
	assert (G.exists[u]);
	assert (v >= 0);
	assert (v < G.z);
	assert (G.exists[v]);

	if (u == v) { // self-loop case
		G.adja[u].push_back(u);
		G.deg[u] += 1;
		if (w == Weighted::weighted) {
			G.edgeWeights[u].push_back(ew);
		}
		for (index attrId = 0; attrId < G.edgeMaps_double.size(); ++attrId) {
			double defaultAttr = G.edgeAttrDefaults_double[attrId];
			G.edgeMaps_double[attrId][u].push_back(defaultAttr);
		}
	} else {
		// set adjacency
		G.adja[u].push_back(v);
		G.adja[v].push_back(u);
		// increment degree counters
		G.deg[u]++;
		G.deg[v]++;
		// set edge weight
		if (w == Weighted::weighted) {
			G.edgeWeights[u].push_back(ew);
			G.edgeWeights[v].push_back(ew);	
		}
		// loop over all attributes, setting default attr
		for (index attrId = 0; attrId < G.edgeMaps_double.size(); ++attrId) {
			double defaultAttr = G.edgeAttrDefaults_double[attrId];
			G.edgeMaps_double[attrId][u].push_back(defaultAttr);
			G.edgeMaps_double[attrId][v].push_back(defaultAttr);
		}
	}

	G.m++; // increasing the number of edges
}

template<Weighted w>
void addEdge_impl(BasicGraph<w, Directed::directed>& G, node u, node v, edgeweight ew) {
	assert (u >= 0);
	assert (u < G.z);
	assert (G.exists[u]);
	assert (v >= 0);
	assert (v < G.z);
	assert (G.exists[v]);

	G.m++; // increasing the number of edges
	G.inDeg[u]++;
	G.outDeg[v]++;

	G.outEdges[u].push_back(v);
	if (w == Weighted::weighted) {
		G.edgeWeights[u].push_back(ew);
	}
	G.inEdges[v].push_back(u);

	// loop over all attributes, setting default attr
	for (index attrId = 0; attrId < G.edgeMaps_double.size(); ++attrId) {
		double defaultAttr = G.edgeAttrDefaults_double[attrId];
		G.edgeMaps_double[attrId][u].push_back(defaultAttr);
		G.edgeMaps_double[attrId][v].push_back(defaultAttr);
	}

}

template<Weighted w>
void removeEdge_impl(const BasicGraph<w, Directed::undirected>& G, node u, node v) {
	index ui = G.find(v, u);
	index vi = G.find(u, v);

	if (vi == none) {
		std::stringstream strm;
		strm << "edge (" << u << "," << v << ") does not exist";
		throw std::runtime_error(strm.str());
		// TODO: what if edge does not exist?
	} else {
		G.adja[u][vi] = none;
		G.adja[v][ui] = none;
		// decrement degree counters
		G.deg[u]--;
		if (u != v) { // self-loops are counted only once
			G.deg[v]--;
		}
		if (w == Weighted::weighted) {
			// remove edge weight
			G.weights[u][vi] = G.nullWeight;
			G.weights[v][ui] = G.nullWeight;
		}

		// TODO: remove attributes

		G.m--; // decreasing the number of edges

	}
}

template<Weighted w>
void removeEdge_impl(const BasicGraph<w, Directed::directed>& G, node u, node v) {
	// remove adjacency
	index ui = G.find(v, u);
	index vi = G.find(u, v);

	if (vi == none) {
		std::stringstream strm;
		strm << "edge (" << u << "," << v << ") does not exist";
		throw std::runtime_error(strm.str());
	} else {
		G.outEdges[u][vi] = none;
		G.inEdges[v][ui] = none;
		// decrement degree counters
		G.deg[u].out--;
		G.deg[v].in--;
		if (w == Weighted::weighted) {
			// remove edge weight
			G.edgeWeights[u][vi] = G.nullWeight;
			G.edgeWeights[v][ui] = G.nullWeight;
		}

		// TODO: remove attributes

		G.m--; // decreasing the number of edges
	}
}

template<Weighted w, Directed d>
bool BasicGraph<w, d>::hasEdge(node u, node v) const {
	return find(u, v) != none;
}

template<Weighted w, Directed d>
node BasicGraph<w, d>::mergeEdge(node u, node v, bool discardSelfLoop) {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** GLOBAL PROPERTIES **/

template<Weighted w, Directed d>
count BasicGraph<w, d>::numberOfSelfLoops() const {
	count c = 0;
	forEdges([&](node u, node v) {
		if (u == v) {
			c += 1;
		}
	});
	return c;
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

template<Directed d>
edgeweight weight_impl(const BasicGraph<Weighted::unweighted, d>& G, node u, node v) {
	return G.hasEdge(u, v) ? defaultEdgeWeight : nullWeight;
}

template<Directed d>
edgeweight weight_impl(const BasicGraph<Weighted::weighted, d>& G, node u, node v) {
	index vi = G.find(u, v);
	if (vi != none) {
		return G.edgeWeights[u][vi];
	} else {
		return nullWeight;
	}
}

template<>
void BasicGraph<Weighted::unweighted, Directed::undirected>::setWeight(node u, node v, edgeweight ew) = delete;

template<>
void BasicGraph<Weighted::weighted, Directed::undirected>::setWeight(node u, node v, edgeweight ew) {
	if (u == v) {
		// self-loop case
		index ui = find(u, u);
		if (ui != none) {
			edgeWeights[u][ui] = ew;
		} else {
			addEdge(u, u, ew);
		}
	} else {
		index vi = find(u, v);
		index ui = find(v, u);
		if ((vi != none) && (ui != none)) {
			edgeWeights[u][vi] = ew;
			edgeWeights[v][ui] = ew;
		} else {
			addEdge(u, v, ew);
		}
	}
}

template<>
void BasicGraph<Weighted::unweighted, Directed::directed>::setWeight(node u, node v, edgeweight ew) = delete;

template<>
void BasicGraph<Weighted::weighted, Directed::directed>::setWeight(node u, node v, edgeweight ew) {
	index vi = find(u, v);
	index ui = find(v, u); // ToDo use findIn
	if ((vi != none) && (ui != none)) {
		edgeWeights[u][vi] = ew;
		edgeWeights[v][ui] = ew;
	} else {
		addEdge(u, v, ew);
	}
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::increaseWeight(node u, node v, edgeweight ew) {
	setWeight(u, v, weight(u, v) + ew);
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
