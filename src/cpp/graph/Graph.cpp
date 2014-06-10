/*
 * BasicGraph.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#include <sstream>

#include "Graph.h"
#include "../auxiliary/Log.h"
 
namespace NetworKit {

/** CONSTRUCTORS **/

Graph::Graph(count n, bool weighted, bool directed) :
	n(n),
	m(0),
	z(n),
	t(0),

	/*shows wether graph is weighted*/
	weighted(weighted), 

	/*shows wether graph is directed*/
	directed(directed),
	exists(n, true),

	/* for directed graphs inDeg stores the incoming degree of a node, for undirected graphs inDeg is not used*/
	inDeg(directed ? n : 0, 0),

	/* for directed graphs outDeg stores the outgoing degree of a node, for undirected graphs outEdges stores the incoming degree of a node*/
	outDeg(n, 0),

	/* for directed graphs inEdges stores an adjacencylist only considering incoming edges, for undirected graphs inEdges is not used*/
	inEdges(directed ? n : 0),

	/* for directed graphs outEdges stores an adjacencylist only considering outgoing edges, for undirected graphs outEdges stores the adjacencylist of 
	undirected edges*/
	outEdges(n), 
	inEdgeWeights(weighted && directed ? n : 0), 
	outEdgeWeights(weighted ? n : 0) {

	// set name from global id
	static count nextGraphId = 1;
	id = nextGraphId++;
	std::stringstream sstm;
	sstm << "G#" << id;
	name = sstm.str();
}

//only to be used by Cython
void Graph::stealFrom(Graph& input) {
	*this = std::move(input);
}


/** PRIVATE HELPERS **/

index Graph::indexInInEdgeArray(node u, node v) const {
	if (!directed) {
		return indexInOutEdgeArray(u, v);
	}
	for (index i = 0; i < inEdges[u].size(); i++) {
		node x = inEdges[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

index Graph::indexInOutEdgeArray(node u, node v) const {
	for (index i = 0; i < outEdges[u].size(); i++) {
		node x = outEdges[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

edgeweight Graph::outEdgeWeightFromIndex(node u, index vi) const {
	if (vi == none) {
		return nullWeight;
	} else {
		return weighted ? outEdgeWeights[u][vi] : defaultEdgeWeight;
	}
}


/** GRAPH INFORMATION **/

std::string Graph::typ() const { 
	if (weighted) {
		return directed ? "WeightedDirectedGraph" : "WeightedGraph";
	} else {
		return directed ? "DirectedGraph" : "Graph";
	}
}

count Graph::getMemoryUsage() const {
	count mem = 0;

	mem += exists.capacity() / 8;

	for (auto& w : inEdgeWeights) {
		mem += sizeof(edgeweight) * w.capacity();
	}
	for (auto& w : outEdgeWeights) {
		mem += sizeof(edgeweight) * w.capacity();
	}

	mem += sizeof(count) * inDeg.capacity();
	mem += sizeof(count) * outDeg.capacity();
	
	for (auto& a : inEdges) {
		mem += sizeof(node) * a.capacity();
	}
	for (auto& a : outEdges) {
		mem += sizeof(node) * a.capacity();
	}

	for (auto& map : edgeMaps_double) {
		for (auto& a : map) {
			mem += sizeof(double) * a.capacity();
		}
	}

	return mem;
}

void Graph::shrinkToFit() {
	exists.shrink_to_fit();

	inEdgeWeights.shrink_to_fit();
	for (auto& w : inEdgeWeights) {
		w.shrink_to_fit();
	}

	outEdgeWeights.shrink_to_fit();
	for (auto& w : outEdgeWeights) {
		w.shrink_to_fit();
	}

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

	edgeMaps_double.shrink_to_fit();
	for (auto& map : edgeMaps_double) {
		map.shrink_to_fit();
		for (auto& a : map) {
			a.shrink_to_fit();
		}
	}

}

std::string Graph::toString() const {
	std::stringstream strm;
	strm << typ() << "(name=" << getName() << ", n=" << numberOfNodes()
			<< ", m=" << numberOfEdges() << ")";
	return strm.str();
}


/** NODE MODIFIERS **/

node Graph::addNode() {
	node v = z;	// node gets maximum id
	z++;	// increment node range
	n++;	// increment node count

	// update per node data structures
	exists.push_back(true);

	// update per node data structures
	if (weighted) {

		std::vector<edgeweight> edgeWeight;
		inEdgeWeights.push_back(edgeWeight);
		outEdgeWeights.push_back(edgeWeight);
	}

	outDeg.push_back(0);
	outEdges.push_back(std::vector<node>{});
	if (directed) {
		inDeg.push_back(0);
		inEdges.push_back(std::vector<node>{});
	}

	// update edge attribute data structures
	for (size_t attrId = 0; attrId < this->edgeMaps_double.size(); attrId++) {
		std::vector<double> attrVector;
		this->edgeMaps_double[attrId].push_back(attrVector);
	}

	return v;
}

node Graph::addNode(float x, float y) {
	node v = addNode();
	std::vector<float> coords = {x, y};
	coordinates.addCoordinates(coords);
	return v;
}

void Graph::removeNode(node v) {
	assert (v < z);
	assert (exists[v]);

	if (!isIsolated(v)) {
		throw std::runtime_error("nodes must be isolated (degree 0) before they can be removed");
	}

	exists[v] = false;
	n--;
}


/** NODE PROPERTIES **/

edgeweight Graph::weightedDegree(node v) const {
	if (weighted) {
		edgeweight sum = 0.0;
		forWeightedNeighborsOf(v, [&](node u, edgeweight ew) {
			sum += ew;
		});
		return sum;
	}
	return defaultEdgeWeight * degree(v);
}

node Graph::randomNode() const {
	if (numberOfNodes() == 0) {
		return none;
	}

	node v;
	do {
		v = Aux::Random::integer(z - 1);
	} while (!exists[v]);

	return v;
}

node Graph::randomNeighbor(node u) const {
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

void Graph::addEdge(node u, node v, edgeweight ew) {
	assert (u >= 0);
	assert (u < z);
	assert (exists[u]);
	assert (v >= 0);
	assert (v < z);
	assert (exists[v]);

	m++; // increase number of edges
	outDeg[u]++;
	outEdges[u].push_back(v);

	if (directed) {
		inDeg[v]++;
		inEdges[v].push_back(u);

		if (weighted) {
			inEdgeWeights[v].push_back(ew);
			outEdgeWeights[u].push_back(ew);
		}

		// loop over all attributes, setting default attr
		for (index attrId = 0; attrId < edgeMaps_double.size(); ++attrId) {
			double defaultAttr = edgeAttrDefaults_double[attrId];
			edgeMaps_double[attrId][u].push_back(defaultAttr);
			edgeMaps_double[attrId][v].push_back(defaultAttr);
		}
	} else if (u == v) { // self-loop case
		if (weighted) {
			outEdgeWeights[u].push_back(ew);
		}

		for (index attrId = 0; attrId < edgeMaps_double.size(); ++attrId) {
			double defaultAttr = edgeAttrDefaults_double[attrId];
			edgeMaps_double[attrId][u].push_back(defaultAttr);
		}
	} else { // undirected, no self-loop
		outDeg[v]++;
		outEdges[v].push_back(u);

		if (weighted) {
			outEdgeWeights[u].push_back(ew);
			outEdgeWeights[v].push_back(ew);
		}

		// loop over all attributes, setting default attr
		for (index attrId = 0; attrId < edgeMaps_double.size(); ++attrId) {
			double defaultAttr = edgeAttrDefaults_double[attrId];
			edgeMaps_double[attrId][u].push_back(defaultAttr);
			edgeMaps_double[attrId][v].push_back(defaultAttr);
		}
	}
}

void Graph::removeEdge(node u, node v) {
	index vi = indexInOutEdgeArray(u, v);
	index ui = indexInInEdgeArray(v, u);

	if (vi == none) {
		std::stringstream strm;
		strm << "edge (" << u << "," << v << ") does not exist";
		throw std::runtime_error(strm.str());
	}

	m--; // decrease number of edges
	outDeg[u]--;
	outEdges[u][vi] = none;
	if (weighted) {
		outEdgeWeights[u][vi] = nullWeight;
	}

	if (directed) {
		inDeg[v]--;
		inEdges[v][ui] = none;
		if (weighted) {
			inEdgeWeights[v][ui] = nullWeight;
		}
	} else if (u != v) {
		// undirected, not self-loop
		outDeg[v]--;
		outEdges[v][ui] = none;
		if (weighted) {
			outEdgeWeights[v][ui] = nullWeight;
		}
	}

	// dose not make a lot of sense do remove attributes,
	// cause the edge is marked as deleted and we have no null values for the attributes
}

bool Graph::hasEdge(node u, node v) const {
	return indexInOutEdgeArray(u, v) != none;
}

std::pair<node, node> Graph::randomEdge() const {
	// TODO this is relativly fast, but not a uniform random edge!
	node u;
	do {
		u = randomNode();
	} while (outDeg[u] == 0);
	node v = randomNeighbor(u);
	return std::make_pair(u, v);
}


/** GLOBAL PROPERTIES **/

count Graph::numberOfSelfLoops() const {
	count c = 0;
	forEdges([&](node u, node v) {
		if (u == v) {
			c += 1;
		}
	});
	return c;
}


/** EDGE ATTRIBUTES **/

edgeweight Graph::weight(node u, node v) const {
	index vi = indexInOutEdgeArray(u, v);
	return outEdgeWeightFromIndex(u, vi);
}

void Graph::setWeight(node u, node v, edgeweight ew) {
	if (!weighted) {
		throw std::runtime_error("Cannot set edge weight in unweighted graph.");
	}

	index vi = indexInOutEdgeArray(u, v);
	if (vi == none) {
		// edge does not exits, create it, but warn user
		WARN("Setting edge weight of an not existing edge will create the edge.");
		addEdge(u, v, ew);
		return;
	}

	outEdgeWeights[u][vi] = ew;
	if (directed) {
		index ui = indexInInEdgeArray(v, u);
		inEdgeWeights[v][ui] = ew;
	} else if (u != v) {
		// undirected and no self-loop
		index ui = indexInOutEdgeArray(u, v);
		outEdgeWeights[v][ui] = ew;
	}
}

void Graph::increaseWeight(node u, node v, edgeweight ew) {
	if (!weighted) {
		throw std::runtime_error("Cannot increase edge weight in unweighted graph.");
	}

	index vi = indexInOutEdgeArray(u, v);
	if (vi == none) {
		// edge does not exits, create it, but warn user
		WARN("Increasing edge weight of an not existing edge will create the edge.");
		addEdge(u, v, ew);
		return;
	}

	outEdgeWeights[u][vi] += ew;
	if (directed) {
		index ui = indexInInEdgeArray(v, u);
		inEdgeWeights[v][ui] += ew;
	} else if (u != v) {
		// undirected and no self-loop
		index ui = indexInOutEdgeArray(u, v);
		outEdgeWeights[v][ui] += ew;
	}	
}

int Graph::addEdgeAttribute_double(double defaultValue) {
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

double Graph::attribute_double(node u, node v, int attrId) const {
	assert (attrId < edgeMaps_double.size());
	index vi = indexInOutEdgeArray(u, v);
	if (vi != none) {
		return edgeMaps_double[attrId][u][vi];
	} else {
		throw std::runtime_error("Edge does not exist. Can't access double attribute.");
	}
}

void Graph::setAttribute_double(node u, node v, int attrId, double attr) {
	// TODO directed?

	if (u == v) {
		// self-loop case
		index ui = indexInOutEdgeArray(u, u);
		if (ui != none) {
			edgeMaps_double.at(attrId)[u][ui] = attr;
		} else {
			throw std::runtime_error("Edge does not exist. Can't set double attribute.");
		}
	} else {
		index vi = indexInOutEdgeArray(u, v);
		index ui = indexInInEdgeArray(v, u);
		if ((vi != none) && (ui != none)) {
			edgeMaps_double[attrId][u][vi] = attr;
			edgeMaps_double[attrId][v][ui] = attr;
		} else {
			throw std::runtime_error("Edge does not exist. Can't set double attribute.");
		}
	}
}


/** SUMS **/

edgeweight Graph::totalEdgeWeight() const {
	if (weighted) {
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

std::vector<node> Graph::nodes() const {
	std::vector<node> nodes;
	nodes.reserve(numberOfNodes());
	this->forNodes([&](node u) {
		nodes.push_back(u);
	});
	return nodes;
}


std::vector<std::pair<node, node> > Graph::edges() const {
	std::vector<std::pair<node, node> > edges;
	edges.reserve(numberOfEdges());
	this->forEdges([&](node u, node v){
		edges.push_back(std::pair<node, node>(u, v));
	});
	return edges;
	
}

std::vector<node> Graph::neighbors(node u) const {
	std::vector<node> neighbors;
	neighbors.reserve(degree(u));
	this->forNeighborsOf(u, [&](node v) {
		neighbors.push_back(v);
	});
	return neighbors;
}

} /* namespace NetworKit */

