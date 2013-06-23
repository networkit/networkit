/*
 *
 *  Created on: 04.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu), Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#include "Graph.h"

namespace NetworKit {

Graph::Graph(): n(0), m(0), z(n), t(0), weighted(false), deg(z, 0), exists(z, true), adja(z), eweights(z) {
	// set name from global id
	static int64_t graphId = 1;
	std::stringstream sstm;
	sstm << "G#" << graphId++;
	this->name = sstm.str();

}

// TODO: z should probably be n-1, but it breaks some tests
Graph::Graph(count n) : n(n), m(0), z(n), t(0), weighted(false), deg(z, 0), exists(z, true), adja(z), eweights(z) {
	// set name from global id
	static int64_t graphId = 1;
	std::stringstream sstm;
	sstm << "G#" << graphId++;
	this->name = sstm.str();
}

Graph::~Graph() {

}


index Graph::find(node u, node v) const {
	for (index vi = 0; vi < this->adja[u].size(); ++vi) {
		node x = this->adja[u][vi];
		if (x == v) {
			return vi;
		}
	}

	return none;
}

void Graph::addEdge(node u, node v, edgeweight weight) {
	assert (u >= 0);
	assert (u <= this->z);
	assert (v >= 0);
	assert (v <= this->z); // node ids must be in range

	if (u == v) { // self-loop case
		this->adja[u].push_back(u);
		this->deg[u] += 1;
		this->eweights[u].push_back(weight);
		for (index attrId = 0; attrId < this->edgeMaps_double.size(); ++attrId) {
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
		for (index attrId = 0; attrId < this->edgeMaps_double.size(); ++attrId) {
			double defaultAttr = this->edgeAttrDefaults_double[attrId];
			this->edgeMaps_double[attrId][u].push_back(defaultAttr);
			this->edgeMaps_double[attrId][v].push_back(defaultAttr);
		}
	}

	m++; // increasing the number of edges
}

void Graph::removeEdge(node u, node v) {
	// remove adjacency
	index ui = find(v, u);
	index vi = find(u, v);

	if (vi == none) {
		throw std::runtime_error("edge does not exist");
		// TODO: what if edge does not exist?
	} else {
		this->adja[u][vi] = none;
		this->adja[v][ui] = none;
		// decrement degree counters
		this->deg[u] -= 1;
		if (u != v) { // self-loops are counted only once
			this->deg[v] -= 1;
		}
		// remove edge weight
		this->eweights[u][vi] = this->nullWeight;
		this->eweights[v][ui] = this->nullWeight;
		// TODO: remove attributes

		m--; // decreasing the number of edges

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
			addEdge(u, u, w);
		}
	} else {
		index vi = find(u, v);
		index ui = find(v, u);
		if ((vi != none) && (ui != none)) {
			this->eweights[u][vi] = w;
			this->eweights[v][ui] = w;
		} else {
			addEdge(u, v, w);
		}
	}

}

bool Graph::hasEdge(node u, node v) const {
	return (find(u, v) != none);
}

node Graph::addNode() {
	node v = this->z;	// node gets maximum id
	this->z += 1;	// increment node range
	this->n += 1;	// increment node count

	//update per node data structures
	this->deg.push_back(0);
	this->exists.push_back(true);

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


node Graph::addNode(float x, float y) {
	node u = addNode();
	std::vector<float> coords = {x, y};
	coordinates.addCoordinates(coords);
	return u;
}


bool Graph::isEmpty() {
	return (n == 0);
}

count Graph::numberOfNodes() const {
	return this->n;
}

count Graph::degree(node v) const {
	assert (v >= 0);
	assert (v <= this->z); // node ids must be in range
	return deg[v];
}

edgeweight Graph::weightedDegree(node v) const {
	// weighted degree as sum over incident edge weight - self-loops are counted once
	edgeweight wDeg = 0.0;
	for (edgeweight w : this->eweights[v]) {
		wDeg += w;
	}
	return wDeg;
}


count Graph::numberOfEdges() const {
	// returns a field which is updated on addEdge() and removeEdge()
	return m;
}

edgeweight Graph::totalEdgeWeight() const {
	edgeweight sum = 0.0;
	this->forWeightedEdges([&](node u, node v, edgeweight ew) {
		sum += ew;
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

count Graph::numberOfSelfLoops() const {
	count nl = 0;
	this->forEdges([&](node u, node v) {
		if (u == v) {
			nl += 1;
		}
	});
	return nl;
}

void Graph::removeNode(node u) {
	if (this->degree(u) > 0) {
		throw std::runtime_error("nodes must have degree 0 before they can be removed");
	}

	this->exists[u] = false;
	this->n -= 1;
}

bool Graph::hasNode(node u) const {
	return this->exists[u];	// exists array determines whether node is present
}


int Graph::addEdgeAttribute_double(double defaultValue) {
	int attrId = this->edgeMaps_double.size();
	std::vector<std::vector<double> > edgeMap;
	edgeMap.resize(this->n);	// create empty vector<attr> for each node
	this->edgeAttrDefaults_double.push_back(defaultValue);

	if (this->numberOfEdges() > 0) {
		forNodes([&] (node v) {
			// create edgeMaps and fill them with default value
			edgeMap[v].resize(adja[v].size());
			fill(edgeMap[v].begin(), edgeMap[v].end(), defaultValue);
		});
//		throw std::runtime_error("TODO: set attributes for already existing edges");
	}
	this->edgeMaps_double.push_back(edgeMap);

	return attrId;
}

std::vector<node> Graph::nodes() {
	std::vector<node> nodes;
	this->forNodes([&](node u){
		nodes.push_back(u);
	});
	return nodes;
}

void Graph::markAsWeighted() {
	this->weighted = true;
}

bool Graph::isMarkedAsWeighted() {
	return this->weighted;
}

std::vector<std::pair<node, node> > Graph::edges() {
	std::vector<std::pair<node, node> > edges;
	this->forEdges([&](node u, node v){
		edges.push_back(std::pair<node, node>(u, v));
	});
	return edges;
}

void Graph::timeStep() {
	this->t += 1;
}

count Graph::time() {
	return this->t;
}


index Graph::upperNodeIdBound() const {
	return this->z;
}


} /* namespace NetworKit */

