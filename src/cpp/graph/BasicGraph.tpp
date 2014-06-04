/*
 * BasicGraph.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#include <stdexcept>
#include <sstream>

#include "BasicGraph.h"
#include "../auxiliary/Random.h"
#include <stack>
#include <queue>

namespace NetworKit {

namespace graph_impl {


/** CONSTRUCTORS **/

template<Weighted w, Directed d>
BasicGraph<w, d>::BasicGraph(count n, bool dummy) :
	WData(n),
	DData(n),
	n(n),
	m(0),
	z(n),
	t(0),
	exists(n, true) {

	if (dummy && w == Weighted::unweighted) {
		throw std::runtime_error("weighted parameter for constructor is deprecated. Graph_T and DirectedGraph_T are always unweighted. Use WeightedGraph_T or WeightedDirectedGraph_T instead.");
	}
	
	// set name from global id
	static count nextGraphId = 1;
	id = nextGraphId++;
	std::stringstream sstm;
	sstm << "G#" << id;
	name = sstm.str();
}

 
//only to be used by Cython
template<Weighted w, Directed d>
void BasicGraph<w, d>::stealFrom(BasicGraph<w, d>& input) {
	*this = std::move(input);
}

/** PRIVATE HELPERS **/

template<>
inline const std::vector<node>& BasicGraph<Weighted::unweighted, Directed::undirected>::adjaIn(node u) const { return adja[u]; }
template<>
inline const std::vector<node>& BasicGraph<Weighted::weighted, Directed::undirected>::adjaIn(node u) const { return adja[u]; }
template<>
inline const std::vector<node>& BasicGraph<Weighted::unweighted, Directed::directed>::adjaIn(node u) const { return inEdges[u]; }
template<>
inline const std::vector<node>& BasicGraph<Weighted::weighted, Directed::directed>::adjaIn(node u) const { return inEdges[u]; }

template<>
inline const std::vector<node>& BasicGraph<Weighted::unweighted, Directed::undirected>::adjaOut(node u) const { return adja[u]; }
template<>
inline const std::vector<node>& BasicGraph<Weighted::weighted, Directed::undirected>::adjaOut(node u) const { return adja[u]; }
template<>
inline const std::vector<node>& BasicGraph<Weighted::unweighted, Directed::directed>::adjaOut(node u) const { return outEdges[u]; }
template<>
inline const std::vector<node>& BasicGraph<Weighted::weighted, Directed::directed>::adjaOut(node u) const { return outEdges[u]; }


/** GRAPH INFORMATION **/

template<>
inline std::string BasicGraph<Weighted::unweighted, Directed::undirected>::typ() const { return "Graph"; }

template<>
inline std::string BasicGraph<Weighted::weighted, Directed::undirected>::typ() const { return "WeightedGraph"; }

template<>
inline std::string BasicGraph<Weighted::unweighted, Directed::directed>::typ() const { return "DirectedGraph"; }

template<>
inline std::string BasicGraph<Weighted::weighted, Directed::directed>::typ() const { return "WeightedDirectedGraph"; }

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
	node v = z;	// node gets maximum id
	z++;	// increment node range
	n++;	// increment node count

	// update per node data structures
	exists.push_back(true);

	// update per node data structures
	WData::addNode();
	DData::addNode();

	// update edge attribute data structures
	for (size_t attrId = 0; attrId < this->edgeMaps_double.size(); attrId++) {
		std::vector<double> attrVector;
		this->edgeMaps_double[attrId].push_back(attrVector);
	}

	return v;
}

template<Weighted w, Directed d>
node BasicGraph<w, d>::addNode(float x, float y) {
	node v = addNode();
	std::vector<float> coords = {x, y};
	coordinates.addCoordinates(coords);
	return v;
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::removeNode(node v) {
	assert (v < z);
	assert (exists[v]);

	if (!isIsolated(v)) {
		throw std::runtime_error("nodes must be isolated (degree 0) before they can be removed");
	}

	exists[v] = false;
	n--;
}


/** NODE PROPERTIES **/

template<Weighted w>
count degreeIn_impl(const BasicGraph<w, Directed::undirected>& G, node v) {
	return G.deg[v];
}

template<Weighted w>
count degreeIn_impl(const BasicGraph<w, Directed::directed>& G, node v) {
	return G.inDeg[v];
}

template<Weighted w>
count degreeOut_impl(const BasicGraph<w, Directed::undirected>& G, node v) {
	return G.deg[v];
}

template<Weighted w>
count degreeOut_impl(const BasicGraph<w, Directed::directed>& G, node v) {
	return G.outDeg[v];
}

template<Weighted w, Directed d>
edgeweight BasicGraph<w, d>::weightedDegree(node v) const {
	if (w == Weighted::weighted) {
		edgeweight sum = 0.0;
		forWeightedNeighborsOf(v, [&](node u, edgeweight ew) {
			sum += ew;
		});
		return sum;
	}
	return defaultEdgeWeight * degree(v);
}

template<Weighted w, Directed d>
node BasicGraph<w, d>::randomNode() const {
	if (numberOfNodes() == 0) {
		return none;
	}

	node v;
	do {
		v = Aux::Random::integer(z - 1);
	} while (!exists[v]);

	return v;
}

template<Weighted w, Directed d>
node BasicGraph<w, d>::randomNeighbor(node u) const {
	if (degree(u) == 0) {
		return none;
	}

	auto neighbors = adjaOut(u);
	node v;
	do {
		v = Aux::Random::integer(neighbors.size() - 1);
	} while (v == none);

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

		G.addEdgeWeight(u, ew);
		
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
		G.addEdgeWeight(u, ew);
		G.addEdgeWeight(v, ew);

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
	G.outDeg[u]++;
	G.inDeg[v]++;

	G.outEdges[u].push_back(v);
	G.inEdges[v].push_back(u);

	G.addEdgeWeight(u, ew);

	// loop over all attributes, setting default attr
	for (index attrId = 0; attrId < G.edgeMaps_double.size(); ++attrId) {
		double defaultAttr = G.edgeAttrDefaults_double[attrId];
		G.edgeMaps_double[attrId][u].push_back(defaultAttr);
		G.edgeMaps_double[attrId][v].push_back(defaultAttr);
	}

}

template<Weighted w>
void removeEdge_impl(BasicGraph<w, Directed::undirected>& G, node u, node v) {
	index ui = G.indexInEdgeArray(v, u);
	index vi = G.indexInEdgeArray(u, v);

	if (vi == none) {
		std::stringstream strm;
		strm << "edge (" << u << "," << v << ") does not exist";
		throw std::runtime_error(strm.str());
	} else {
		G.m--; // decreasing the number of edges
		G.adja[u][vi] = none;
		G.adja[v][ui] = none;
		// decrement degree counters
		G.deg[u]--;
		if (u != v) { // self-loops are counted only once
			G.deg[v]--;
		}

		G.setEdgeWeight(u, vi, nullWeight);
		G.setEdgeWeight(v, ui, nullWeight);

		// dose not make a lot of sense do remove attributes,
		// cause the edge is marked as deleted and we have no null values for the attributes
	}
}

template<Weighted w>
void removeEdge_impl(BasicGraph<w, Directed::directed>& G, node u, node v) {
	// remove adjacency, for u it is on outgoing edge, at v an incoming
	index ui = G.indexInInEdgeArray(v, u);
	index vi = G.indexInOutEdgeArray(u, v);

	if (vi == none) {
		std::stringstream strm;
		strm << "edge (" << u << "," << v << ") does not exist";
		throw std::runtime_error(strm.str());
	} else {
		G.m--; // decreasing the number of edges
		G.outEdges[u][vi] = none;
		G.inEdges[v][ui] = none;
		// decrement degree counters
		G.inDeg[u]--;
		G.outDeg[v]--;

		G.setEdgeWeight(u, vi, nullWeight);

		// dose not make a lot of sense do remove attributes,
		// cause the edge is marked as deleted and we have no null values for the attributes
	}
}

template<Weighted w, Directed d>
bool BasicGraph<w, d>::hasEdge(node u, node v) const {
	return DData::indexInEdgeArray(u, v) != none;
}

template<Weighted w, Directed d>
node BasicGraph<w, d>::mergeEdge(node u, node v, bool discardSelfLoop) {
	throw std::runtime_error("TODO: implement mergeEdge");
	// if (u != v) {
	// 	node newNode = addNode();

	// 	// self-loop if necessary
	// 	if (! discardSelfLoop) {
	// 		edgeweight selfLoopWeight = this->weight(u, u) + this->weight(v, v) + this->weight(u, v);
	// 		this->addEdge(newNode, newNode, selfLoopWeight);
	// 	}

	// 	// rewire edges from u to newNode
	// 	this->forWeightedEdgesOf(u, [&](node u, node neighbor, edgeweight weight) {
	// 		if (neighbor != u && neighbor != v) {
	// 			this->increaseWeight(neighbor, newNode, this->weight(u, neighbor)); // TODO: make faster
	// 		}
	// 	});

	// 	// rewire edges from v to newNode
	// 	this->forWeightedEdgesOf(v, [&](node v, node neighbor, edgeweight weight) {
	// 		if (neighbor != v && neighbor != u) {
	// 			this->increaseWeight(neighbor, newNode, this->weight(v, neighbor));  // TODO: make faster
	// 		}
	// 	});

	// 	// delete edges of nodes to delete
	// 	this->forEdgesOf(u, [&](node u, node neighbor) {
	// 		this->removeEdge(u, neighbor);
	// 	});
	// 	this->forEdgesOf(v, [&](node v, node neighbor) {
	// 		this->removeEdge(v, neighbor);
	// 	});


	// 	// delete nodes
	// 	this->removeNode(u);
	// 	this->removeNode(v);


	// 	return newNode;
	// }

	// no new node created
	return none;
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


/** EDGE ATTRIBUTES **/

template<>
inline edgeweight BasicGraph<Weighted::unweighted, Directed::undirected>::weight(node u, node v) const { return hasEdge(u, v) ? defaultEdgeWeight : nullWeight; }
template<>
inline edgeweight BasicGraph<Weighted::weighted, Directed::undirected>::weight(node u, node v) const { return edgeWeightFromIndex(u, indexInEdgeArray(u, v)); }
template<>
inline edgeweight BasicGraph<Weighted::unweighted, Directed::directed>::weight(node u, node v) const { return hasEdge(u, v) ? defaultEdgeWeight : nullWeight; }
template<>
inline edgeweight BasicGraph<Weighted::weighted, Directed::directed>::weight(node u, node v) const { return edgeWeightFromIndex(u, indexInEdgeArray(u, v)); }

template<>
inline void BasicGraph<Weighted::unweighted, Directed::undirected>::setWeight(node u, node v, edgeweight ew) {}; // not marked 'deleted' for better compatibility

template<>
inline void BasicGraph<Weighted::weighted, Directed::undirected>::setWeight(node u, node v, edgeweight ew) {
	if (u == v) {
		// self-loop case
		index ui = indexInEdgeArray(u, u);
		if (ui != none) {
			edgeWeights[u][ui] = ew;
		} else {
			addEdge(u, u, ew);
		}
	} else {
		index vi = indexInEdgeArray(u, v);
		index ui = indexInEdgeArray(v, u);
		if ((vi != none) && (ui != none)) {
			edgeWeights[u][vi] = ew;
			edgeWeights[v][ui] = ew;
		} else {
			addEdge(u, v, ew);
		}
	}
}

template<>
inline void BasicGraph<Weighted::unweighted, Directed::directed>::setWeight(node u, node v, edgeweight ew) {}; // not marked 'deleted' for better compatibility

template<>
inline void BasicGraph<Weighted::weighted, Directed::directed>::setWeight(node u, node v, edgeweight ew) {
	index vi = indexInOutEdgeArray(u, v);
	index ui = indexInInEdgeArray(v, u);
	if ((vi != none) && (ui != none)) {
		edgeWeights[u][vi] = ew;
		edgeWeights[v][ui] = ew;
	} else {
		addEdge(u, v, ew);
	}
}

template<>
inline void BasicGraph<Weighted::unweighted, Directed::undirected>::increaseWeight(node u, node v, edgeweight ew) {}; // not marked 'deleted' for better compatibility
template<>
inline void BasicGraph<Weighted::weighted, Directed::undirected>::increaseWeight(node u, node v, edgeweight ew) { setWeight(u, v, weight(u, v) + ew); }
template<>
inline void BasicGraph<Weighted::unweighted, Directed::directed>::increaseWeight(node u, node v, edgeweight ew) {}; // not marked 'deleted' for better compatibility
template<>
inline void BasicGraph<Weighted::weighted, Directed::directed>::increaseWeight(node u, node v, edgeweight ew) { setWeight(u, v, weight(u, v) + ew); }

template<Weighted w, Directed d>
int BasicGraph<w, d>::addEdgeAttribute_double(double defaultValue) {
	int attrId = edgeMaps_double.size();

	std::vector<std::vector<double> > edgeMap(z);
	if (numberOfEdges() > 0) {
		forNodes([&] (node v) {
			// create edgeMaps and fill them with default value
			edgeMap[v].resize(degree(v));
			fill(edgeMap[v].begin(), edgeMap[v].end(), defaultValue);
		});
	}
	
	edgeMaps_double.push_back(edgeMap);
	edgeAttrDefaults_double.push_back(defaultValue);

	return attrId;
}

template<Weighted w, Directed d>
double BasicGraph<w, d>::attribute_double(node u, node v, int attrId) const {
	assert (attrId < edgeMaps_double.size());
	index vi = BasicGraph<w, d>::indexInEdgeArray(u, v);
	if (vi != none) {
		return edgeMaps_double[attrId][u][vi];
	} else {
		throw std::runtime_error("Edge does not exist. Can't access double attribute.");
	}
}

template<Weighted w>
void setAttribute_double_impl(BasicGraph<w, Directed::undirected>& G, node u, node v, int attrId, double attr) {
	if (u == v) {
		// self-loop case
		index ui = G.indexInEdgeArray(u, u);
		if (ui != none) {
			G.edgeMaps_double.at(attrId)[u][ui] = attr;
		} else {
			throw std::runtime_error("Edge does not exist. Can't set double attribute.");
		}
	} else {
		index vi = G.indexInEdgeArray(u, v);
		index ui = G.indexInEdgeArray(v, u);
		if ((vi != none) && (ui != none)) {
			G.edgeMaps_double[attrId][u][vi] = attr;
			G.edgeMaps_double[attrId][v][ui] = attr;
		} else {
			throw std::runtime_error("Edge does not exist. Can't set double attribute.");
		}
	}
}

template<Weighted w>
void setAttribute_double_impl(BasicGraph<w, Directed::directed>& G, node u, node v, int attrId, double attr) {
	index vi = G.indexInOutEdgeArray(u, v);
	index ui = G.indexInInEdgeArray(v, u);
	if ((vi != none) && (ui != none)) {
		G.edgeMaps_double[attrId][u][vi] = attr;
		G.edgeMaps_double[attrId][v][ui] = attr;
	} else {
		throw std::runtime_error("Edge does not exist. Can't set double attribute.");
	}
}


/** SUMS **/

template<Weighted w, Directed d>
edgeweight BasicGraph<w, d>::totalEdgeWeight() const {
	if (Weighted::weighted == w) {
		edgeweight sum = 0.0;
		forWeightedEdges([&](node u, node v, edgeweight ew) {
			sum += ew;
		});
		return sum;
	} else {
		return numberOfEdges() * defaultEdgeWeight;
	}
}


/** Collections **/

template<Weighted w, Directed d>
std::vector<node> BasicGraph<w, d>::nodes() const {
	std::vector<node> nodes;
	nodes.reserve(numberOfNodes());
	this->forNodes([&](node u) {
		nodes.push_back(u);
	});
	return nodes;
}


template<Weighted w, Directed d>
std::vector<std::pair<node, node> > BasicGraph<w, d>::edges() const {
	std::vector<std::pair<node, node> > edges;
	edges.reserve(numberOfEdges());
	this->forEdges([&](node u, node v){
		edges.push_back(std::pair<node, node>(u, v));
	});
	return edges;
	
}


template<Weighted w, Directed d>
std::vector<node> BasicGraph<w, d>::neighbors(node u) const {
	std::vector<node> neighbors;
	neighbors.reserve(degree(u));
	this->forNeighborsOf(u, [&](node v) {
		neighbors.push_back(v);
	});
	return neighbors;
}


/** NODE ITERATORS **/

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forNodes(L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::parallelForNodes(L handle) const {
	#pragma omp parallel for
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<Weighted w, Directed d>
template<typename C, typename L>
void BasicGraph<w, d>::forNodesWhile(C condition, L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			if (!condition()) {
				break;
			}
			handle(v);
		}
	}
}

template<Weighted w, Directed d>
template<typename C, typename L>
void BasicGraph<w, d>::forNodes(C condition, L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			if (condition()) {
				break;
			}
			handle(v);
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forNodesInRandomOrder(L handle) const {
	std::vector<node> randVec = nodes();
	random_shuffle(randVec.begin(), randVec.end());
	for (node v : randVec) {
		handle(v);
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::balancedParallelForNodes(L handle) const {
	#pragma omp parallel for schedule(guided) // TODO: define min block size (and test it!)
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forNodePairs(L handle) const {
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				if (exists[v]) {
					handle(u, v);
				}
			}
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::parallelForNodePairs(L handle) const {
	#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				if (exists[v]) {
					handle(u, v);
				}
			}
		}
	}
}

	
/** EDGE ITERATORS **/

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forEdges(L handle) const {
	for (node u = 0; u < z; ++u) {
		auto& neighbors = adjaOut(u);
		for (node v : neighbors) {
			if (d == Directed::directed) {
				if (v != none) {
					handle(u, v);
				}
			} else {
				// undirected, do not iterate over edges twice
				// {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				if (u >= v) {
					handle(u, v);
				}
			}
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::parallelForEdges(L handle) const {
	#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		auto& neighbors = adjaOut(u);
		for (node v : neighbors) {
			if (d == Directed::directed) {
				if (v != none) {
					handle(u, v);
				}
			} else {
				// undirected, do not iterate over edges twice
				// {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				if (u >= v) {
					handle(u, v);
				}
			}
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forWeightedEdges(L handle) const {
	for (node u = 0; u <z; ++u) {
		auto& neighbors = adjaOut(u);
		for (index i = 0; i < neighbors.size(); ++i) {
			node v = neighbors[i];
			if (d == Directed::directed) {
				if (v != none) {
					edgeweight ew = WData::edgeWeightFromIndex(u, i);
					handle(u, v, ew);
				}
			} else {
				// undirected, do not iterate over edges twice
				// {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				if (u >= v) {
					edgeweight ew = WData::edgeWeightFromIndex(u, i);
					handle(u, v, ew);
				}
			}
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::parallelForWeightedEdges(L handle) const {
	#pragma omp parallel for
	for (node u = 0; u <z; ++u) {
		auto& neighbors = adjaOut(u);
		for (index i = 0; i < neighbors.size(); ++i) {
			node v = neighbors[i];
			if (d == Directed::directed) {
				if (v != none) {
					edgeweight ew = WData::edgeWeightFromIndex(u, i);
					handle(u, v, ew);
				}
			} else {
				// undirected, do not iterate over edges twice
				// {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				if (u >= v) {
					edgeweight ew = WData::edgeWeightFromIndex(u, i);
					handle(u, v, ew);
				}
			}
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forEdgesWithAttribute_double(int attrId, L handle) const {
	for (node u = 0; u < z; ++u) {
		auto& neighbors = adjaOut(u);
		for (index i = 0; i < neighbors.size(); ++i) {
			node v = neighbors[i];
			if (d == Directed::directed) {
				if (v != none) {
					double attr = edgeMaps_double[attrId][u][i];
					handle(u, v, attr);
				}
			} else {
				// undirected, do not iterate over edges twice
				// {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				if (u >= v) {
					double attr = edgeMaps_double[attrId][u][i];
					handle(u, v, attr);
				}
			}
		}
	}
}


/** NEIGHBORHOOD ITERATORS **/

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forNeighborsOf(node u, L handle) const {
	auto& neighbors = adjaOut(u);
	for (auto v : neighbors) {
		if (v != none) {
			handle(v);
		}
	}	
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forWeightedNeighborsOf(node u, L handle) const {
	auto& neighbors = adjaOut(u);
	for (index i = 0; i < neighbors.size(); i++) {
		node v = neighbors[i];
		if (v != none) {
			edgeweight ew = WData::edgeWeightFromIndex(u, i);
			handle(v, ew);
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forEdgesOf(node u, L handle) const {
	auto& neighbors = adjaOut(u);
	for (auto v : neighbors) {
		if (v != none) {
			handle(u, v);
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forWeightedEdgesOf(node u, L handle) const {
	auto& neighbors = adjaOut(u);
	for (index i = 0; i < neighbors.size(); i++) {
		node v = neighbors[i];
		if (v != none) {
			edgeweight ew = WData::edgeWeightFromIndex(u, i);
			handle(u, v, ew);
		}
	}
}


/** REDUCTION ITERATORS **/

template<Weighted w, Directed d>
template<typename L>
double BasicGraph<w, d>::parallelSumForNodes(L handle) const {
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			sum += handle(v);
		}
	}
	return sum;
}

template<Weighted w, typename L>
double parallelSumForNodes_impl(const BasicGraph<w, Directed::undirected>& G, L handle) {
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (node u = 0; u < G.z; u++) {
		for (index i = 0; i < G.adja[u].size(); i++) {
			node v = G.adja[u][i];
			edgeweight ew = (w == Weighted::weighted) ? G.edgeWeights[u][i] : defaultEdgeWeight;
			if (v != none) {
				sum += handle(u, v, ew);
			}
		}
	}
	return sum;
}

template<Weighted w, typename L>
double parallelSumForNodes_impl(const BasicGraph<w, Directed::directed>& G, L handle) {
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (node u = 0; u < G.z; u++) {
		for (index i = 0; i < G.outEdges[u].size(); i++) {
			node v = G.outEdges[u][i];
			edgeweight ew = (w == Weighted::weighted) ? G.edgeWeights[u][i] : defaultEdgeWeight;
			if (v != none) {
				sum += handle(u, v, ew);
			}
		}
	}
	return sum;
}


/** GRAPH SEARCHES **/

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::BFSfrom(node r, L handle) const {
	std::vector<bool> marked(z);
	std::queue<node> q;
	q.push(r); // enqueue root
	marked[r] = true;
	do {
		node u = q.front();
		q.pop();
		// apply function
		handle(u);
		forNeighborsOf(u, [&](node v) {
			if (!marked[v]) {
				q.push(v);
				marked[v] = true;
			}
		});
	} while (!q.empty());
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::BFSEdgesfrom(node r, L handle) const {
	std::vector<bool> marked(z);
	std::queue<node> q;
	q.push(r); // enqueue root
	marked[r] = true;
	do {
		node u = q.front();
		q.pop();
		// apply function
		forNeighborsOf(u, [&](node v) {
			if (!marked[v]) {
				handle(u, v);
				q.push(v);
				marked[v] = true;
			}
		});
	} while (!q.empty());
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::DFSfrom(node r, L handle) const {
	std::vector<bool> marked(z);
	std::stack<node> s;
	s.push(r); // enqueue root
	marked[r] = true;
	do {
		node u = s.top();
		s.pop();
		// apply function
		handle(u);
		forNeighborsOf(u, [&](node v) {
			if (!marked[v]) {
				s.push(v);
				marked[v] = true;
			}
		});
	} while (!s.empty()); 
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::DFSEdgesfrom(node r, L handle) const {
	std::vector<bool> marked(z);
	std::stack<node> s;
	s.push(r); // enqueue root
	marked[r] = true;
	do {
		node u = s.top();
		s.pop();
		// apply function
		forNeighborsOf(u, [&](node v) {
			if (!marked[v]) {
				handle(u, v);
				s.push(v);
				marked[v] = true;
			}
		});
	} while (!s.empty()); 
}

} /* namespace graph_impl */

} /* namespace NetworKit */

