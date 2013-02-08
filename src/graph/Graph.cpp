/*
 * Graph.cpp
 *
 *  Created on: 04.02.2013
 *      Author: cls
 */

#include "Graph.h"

namespace EnsembleClustering {

Graph::Graph(count n) : n(n), deg(n, 0), adja(n), eweights(n) {
	// set name from global id
	static int64_t graphId = 1;
	std::stringstream sstm;
	sstm << "G#" << graphId++;
	this->name = sstm.str();
}

Graph::~Graph() {
	// TODO Auto-generated destructor stub
}


// TODO: replace by for_each
index Graph::find(node u, node v) const {
	index vi = 0;
	for (node x : this->adja[u]) {
		if (x == v) {
			return vi;
		}
		vi++;
	}

	return none;
}

void Graph::insertEdge(node u, node v, edgeweight weight) {
	if (u == v) { // self-loop case
		this->adja[u].push_back(u);
		this->deg[u] += 1;
		this->eweights[u].push_back(weight);
		for (int attrId = 0; attrId < this->edgeMaps_double.size(); ++attrId) {
			double defaultAttr = this->edgeAttrDefaults_double[attrId];
			this->edgeMaps_double[attrId][u].push_back(defaultAttr);
		}
	} else {
		// set adjacency
		this->adja[u].push_back(v);
		this->adja[v].push_back(u);
		// increment degree counters
		this->deg[u] += 1;
		this->deg[v] += 1;
		// set edge weight
		this->eweights[u].push_back(weight);
		this->eweights[v].push_back(weight);
		// loop over all attributes, setting default attr
		for (int attrId = 0; attrId < this->edgeMaps_double.size(); ++attrId) {
			double defaultAttr = this->edgeAttrDefaults_double[attrId];
			this->edgeMaps_double[attrId][u].push_back(defaultAttr);
			this->edgeMaps_double[attrId][v].push_back(defaultAttr);
		}
	}
}

void Graph::removeEdge(node u, node v) {
	// remove adjacency
	index vi = find(u, v);
	index ui = find(v, u);
	if (vi == none) {
		throw std::runtime_error("edge does not exist");
		// TODO: what if edge does not exist?
	} else {
		this->adja[u][vi] = none;
		this->adja[v][ui] = none;
		// decrement degree counters
		this->deg[u] -= 1;
		this->deg[v] -= 1;
		// remove edge weight
		this->eweights[u][vi] = this->nullWeight;
		this->eweights[v][ui] = this->nullWeight;
		// TODO: remove attributes

	}
}

edgeweight Graph::weight(node u, node v) const {
	index vi = find(u, v);
	if (vi != none) {
		return this->eweights[u][vi];
	} else {
		return 0.0;
	}
}

void Graph::setWeight(node u, node v, edgeweight w) {
	if (u == v) { 		// self-loop case
		index ui = find(u, u);
		if (ui != none) {
			this->eweights[u][ui] = w;
		} else {
			insertEdge(u, u, w);
		}
	} else {
		index vi = find(u, v);
		index ui = find(v, u);
		if ((vi != none) && (ui != none)) {
			this->eweights[u][vi] = w;
			this->eweights[v][ui] = w;
		} else {
			insertEdge(u, v, w);
		}
	}

}

bool Graph::hasEdge(node u, node v) const {
	return (find(u, v) != none);
}

node Graph::addNode() {
	node v = this->n;
	this->n += 1;

	//update per node data structures
	this->deg.push_back(0);

	// update per edge data structures
	std::vector<node> adjacencyVector;	// vector of adjacencies for new node
	std::vector<edgeweight> edgeWeightVector;	// vector of edge weights for new node
	this->adja.push_back(adjacencyVector);
	this->eweights.push_back(edgeWeightVector);

	// update edge attribute data structures
	for (int attrId = 0; attrId < this->edgeMaps_double.size(); ++attrId) {
		std::vector<double> attrVector;
		this->edgeMaps_double[attrId].push_back(attrVector);
	}

	return v;
}

void Graph::extendNodeRange(int64_t n) {
	throw std::runtime_error("TODO");
	// TODO:
}

bool Graph::isEmpty() {
	return (n == 0);
}

int64_t Graph::numberOfNodes() const {
	return this->n;
}

count Graph::degree(node v) const {
	return deg[v];
}

edgeweight Graph::weightedDegree(node v) const {
	// weighted degree as sum over incident edge weight
	edgeweight wDeg = 0.0;
	for (edgeweight w : this->eweights[v]) {
		wDeg += w;
	}
	return wDeg;
}



int64_t Graph::numberOfEdges() const {
	// sum over all stored degrees
	// TODO: parallel sum?
	count mm = 0;
	this->forNodes([&](node v) {
		mm += this->deg[v];
	});
	count m = mm / 2;
	return m;
}

edgeweight Graph::totalEdgeWeight() {
	edgeweight sum = 0.0;
	sum = this->parallelSumForWeightedEdges([&](node u, node v, edgeweight ew) {
		return ew;
	});
	return sum;
}

void Graph::setName(std::string name) {
	this->name = name;
}

std::string Graph::toString() {
	std::stringstream strm;
	strm << "Graph(name=" << this->getName() << ", n=" << this->numberOfNodes()
			<< ", m=" << this->numberOfEdges() << ")";
	return strm.str();
}

std::string Graph::getName() {
	// TODO: unneccessary if name becomes public attribute
	return this->name;
}

edgeweight Graph::totalNodeWeight() {
	throw std::runtime_error("DEPRECATED");
}

void Graph::setAttribute_double(node u, node v, int attrId, double attr) {


	if (u == v) { 		// self-loop case
		index ui = find(u, u);
		if (ui != none) {
			this->edgeMaps_double.at(attrId)[u][ui] = attr;
		} else {
			throw std::runtime_error("What if edge does not exist?");
		}
	} else {
		index vi = find(u, v);
		index ui = find(v, u);
		if ((vi != none) && (ui != none)) {
			// DEBUG
			int s = this->edgeMaps_double.size();
			int sm = this->edgeMaps_double[attrId].size();
			int smu = this->edgeMaps_double[attrId][u].size();
			int smv = this->edgeMaps_double[attrId][v].size();
			// DEBUG

			this->edgeMaps_double[attrId][u][vi] = attr;
			this->edgeMaps_double[attrId][v][ui] = attr;
		} else {
			throw std::runtime_error("What if edge does not exist?");
		}
	}
}

double Graph::attribute_double(node u, node v, int attrId) const {
	assert (attrId < this->edgeMaps_double.size());
	index vi = find(u, v);
	if (vi != none) {
		return this->edgeMaps_double[attrId][u][vi];
	} else {
		throw std::runtime_error("TODO: what if edge does not exist?");
	}

}

int Graph::addEdgeAttribute_double(double defaultValue) {
	int attrId = this->edgeMaps_double.size();
	std::vector<std::vector<double> > edgeMap;
	edgeMap.resize(this->n);	// create empty vector<attr> for each node
	this->edgeMaps_double.push_back(edgeMap);
	this->edgeAttrDefaults_double.push_back(defaultValue);

	if (this->numberOfEdges() > 0) {
		throw std::runtime_error("TODO: set attributes for already existing edges");
	}
	// TODO: set attribute for already existing edges
	return attrId;
}

} /* namespace EnsembleClustering */




