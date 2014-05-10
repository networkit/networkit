/*
 * AbstractGraph.cpp
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#include "DirectedGraph.h"

namespace NetworKit {


// TODO: z should probably be n-1, but it breaks some tests
DirectedGraph::DirectedGraph(count n, bool weighted):
	AbstractGraph(n, weighted),
	deg(z),
	adja(z),
	inOut(z, 0) {
	if (weighted) eweights.resize(z); // edge weight array is only initialized for weighted graphs
}

DirectedGraph::~DirectedGraph() {

}

//only to be used by Cython
void DirectedGraph::stealFrom(DirectedGraph& input) {
	*this = std::move(input);
}


index DirectedGraph::findIn(node u, node v) const {
	for (index i = 0; i < this->inOut[u]; i++) {
		node x = this->adja[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

index DirectedGraph::findOut(node u, node v) const {
	for (index i = this->inOut[u]; i < this->adja[u].size(); i++) {
		node x = this->adja[u][i];
		if (x == v) {
			return i;
		}
	}

	return none;
}


/** NODE MODIFIERS **/

node DirectedGraph::addNode() {
	node v = AbstractGraph::addNode();

	//update per node data structures
	this->deg.push_back(NodeDegree{});

	// update per edge data structures
	this->adja.push_back(_Vector<node>{});
	if (weighted) {
		// vector of edge weights for new node
		this->eweights.push_back(_Vector<edgeweight>{});
	}

	// TODO: support attributes in directed graph
	// update edge attribute data structures
	// for (int attrId = 0; attrId < (int) this->edgeMaps_double.size(); ++attrId) {
	// 	std::vector<double> attrVector;
	// 	this->edgeMaps_double[attrId].push_back(attrVector);
	// }

	return v;
}

node DirectedGraph::addNode(float x, float y) {
	node v = addNode();
	std::vector<float> coords = {x, y};
	coordinates.addCoordinates(coords);
	return v;
}


/** NODE PROPERTIES **/

count DirectedGraph::degree(node v) const {
	assert (v >= 0);
	assert (v < this->z); // node ids must be in range
	return this->deg[v].total();
}

count DirectedGraph::degreeIn(node v) const {
	assert (v >= 0);
	assert (v < this->z); // node ids must be in range
	return this->deg[v].in;
}

count DirectedGraph::degreeOut(node v) const {
	assert (v >= 0);
	assert (v < this->z); // node ids must be in range
	return this->deg[v].out;
}


/** EDGE MODIFIERS **/

void DirectedGraph::addEdge(node u, node v, edgeweight weight) {
	assert (u >= 0);
	assert (u < this->z);
	assert (v >= 0);
	assert (v < this->z);

	m++; // increasing the number of edges
	this->deg[u].out++;
	this->deg[v].in++;

	this->adja[u].push_back(v);
	if (weighted) {
		this->eweights[u].push_back(weight);
	}


	if (this->deg[v].out > 0) {
		// move first outgoing edges to the end of the vector and insert incoming edge
		index i = this->inOut[v];
		node w = this->adja[v][i];
		this->adja[v].push_back(w);
		this->adja[v][i] = u;
		if (weighted) {
			edgeweight we = this->eweights[v][i];
			this->eweights[v].push_back(we);
			this->eweights[v][i] = weight;
		}
	} else {
		this->adja[v].push_back(u);
		if (weighted) {
			this->eweights[v].push_back(weight);
		}
	}
	this->inOut[v]++;

	// TODO: support attributes in directed graph
	// // loop over all attributes, setting default attr
	// for (index attrId = 0; attrId < this->edgeMaps_double.size(); ++attrId) {
	// 	double defaultAttr = this->edgeAttrDefaults_double[attrId];
	// 	this->edgeMaps_double[attrId][u].push_back(defaultAttr);
	// 	this->edgeMaps_double[attrId][v].push_back(defaultAttr);
	// }
}

// void DirectedGraph::removeEdge(node u, node v) {
// 	// remove adjacency
// 	index ui = find(v, u);
// 	index vi = find(u, v);

// 	if (vi == none) {
// 		std::stringstream strm;
// 		strm << "edge (" << u << "," << v << ") does not exist";
// 		throw std::runtime_error(strm.str());
// 		// TODO: what if edge does not exist?
// 	} else {
// 		this->adja[u][vi] = none;
// 		this->adja[v][ui] = none;
// 		// decrement degree counters
// 		this->deg[u] -= 1;
// 		if (u != v) { // self-loops are counted only once
// 			this->deg[v] -= 1;
// 		}
// 		if (weighted) {
// 			// remove edge weight
// 			this->eweights[u][vi] = this->nullWeight;
// 			this->eweights[v][ui] = this->nullWeight;
// 		}

// 		// TODO: remove attributes

// 		m--; // decreasing the number of edges

// 	}
// }


/** GLOBAL PROPERTIES **/

bool DirectedGraph::isDirected() const {
	return true;
}


} /* namespace NetworKit */
