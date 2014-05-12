/*
 *
 *  Created on: 04.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu), Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#include "Graph.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

Graph::Graph(count n, bool weighted):
	AbstractGraph(n, weighted),
	deg(z, 0),
	adja(z) {
	if (weighted) {
		eweights.resize(z); // edge weight array is only initialized for weighted graphs
	}
}

Graph::~Graph() {

}

count Graph::getMemoryUsage() const {
	count mem = AbstractGraph::getMemoryUsage();
	mem += sizeof(deg) + sizeof(count) * deg.capacity();
	
	mem += sizeof(adja) + sizeof(_Vector<node>) * adja.capacity();
	for (auto& a : adja) {
		mem += sizeof(node) * a.capacity();
	}

	mem += sizeof(eweights) + sizeof(_Vector<edgeweight>) * eweights.capacity();
	for (auto& w : eweights) {
		mem += sizeof(edgeweight) * w.capacity();
	}

	// TODO attribute stuff

	return mem;
}

void Graph::shrinkToFit() {
	AbstractGraph::shrinkToFit();
	
	deg.shrink_to_fit();
	adja.shrink_to_fit();
	for (auto& a : adja) {
		a.shrink_to_fit();
	}
	eweights.shrink_to_fit();
	for (auto& w : eweights) {
		w.shrink_to_fit();
	}

	// TODO attribute stuff
}

//only to be used by Cython
void Graph::stealFrom(Graph& input) {
	*this = std::move(input);
}


index Graph::find(node u, node v) const {
	for (index i = 0; i < this->adja[u].size(); i++) {
		node x = this->adja[u][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}


/** NODE MODIFIERS **/

node Graph::addNode() {
	node v = AbstractGraph::addNode();

	//update per node data structures
	this->deg.push_back(0);

	// update per edge data structures
	_Vector<node> adjacencyVector;	// vector of adjacencies for new node
	this->adja.push_back(adjacencyVector);
	if (weighted) {
		_Vector<edgeweight> edgeWeightVector;	// vector of edge weights for new node
		this->eweights.push_back(edgeWeightVector);
	}

	// update edge attribute data structures
	for (int attrId = 0; attrId < (int) this->edgeMaps_double.size(); ++attrId) {
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


/** NODE PROPERTIES **/

bool Graph::isIsolated(node v) const {
	return this->deg[v] == 0;
}

count Graph::degree(node v) const {
	assert (v >= 0);
	assert (v < this->z);
	assert (exists[v]);
	return this->deg[v];
}

count Graph::minDegree() const {
	return *std::min_element(deg.begin(), deg.end());
}

node Graph::argminDegree() const {
	return std::distance(deg.begin(), std::min_element(deg.begin(), deg.end()));
}

count Graph::maxDegree() const {
	return *std::max_element(deg.begin(), deg.end());
}

node Graph::argmaxDegree() const {
	return std::distance(deg.begin(), std::max_element(deg.begin(), deg.end()));
}

edgeweight Graph::weightedDegree(node v) const {
	if (weighted) {
		// weighted degree as sum over incident edge weight - self-loops are counted once
		edgeweight wDeg = 0.0;
		for (edgeweight w : this->eweights[v]) {
			wDeg += w;
		}
		return wDeg;
	} else {
		return this->degree(v);
	}

}

edgeweight Graph::volume(node v) const {
	edgeweight vol = this->weightedDegree(v);
	vol += this->weight(v, v);
	return vol;
}

node Graph::randomNeighbor(node v) const {
	if (degree(v) == 0) {
		TRACE("random neighbor: none");
		return none;
	}

	node randNeigh;
	do {
		randNeigh = adja[v][Aux::Random::integer(adja[v].size())];
	} while (randNeigh == none);

	return randNeigh;
}


/** EDGE MODIFIERS **/

void Graph::addEdge(node u, node v, edgeweight weight) {
	assert (u >= 0);
	assert (u < this->z);
	assert (this->exists[u]);
	assert (v >= 0);
	assert (v < this->z);
	assert (this->exists[v]);

	if (u == v) { // self-loop case
		this->adja[u].push_back(u);
		this->deg[u] += 1;
		if (weighted) this->eweights[u].push_back(weight);
		for (index attrId = 0; attrId < this->edgeMaps_double.size(); ++attrId) {
			double defaultAttr = this->edgeAttrDefaults_double[attrId];
			this->edgeMaps_double[attrId][u].push_back(defaultAttr);
		}
	} else {
		// set adjacency
		this->adja[u].push_back(v);
		this->adja[v].push_back(u);
		// increment degree counters
		this->deg[u]++;
		this->deg[v]++;
		// set edge weight
		if (weighted) {
			this->eweights[u].push_back(weight);
			this->eweights[v].push_back(weight);	
		}
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
		std::stringstream strm;
		strm << "edge (" << u << "," << v << ") does not exist";
		throw std::runtime_error(strm.str());
		// TODO: what if edge does not exist?
	} else {
		this->adja[u][vi] = none;
		this->adja[v][ui] = none;
		// decrement degree counters
		this->deg[u]--;
		if (u != v) { // self-loops are counted only once
			this->deg[v]--;
		}
		if (weighted) {
			// remove edge weight
			this->eweights[u][vi] = this->nullWeight;
			this->eweights[v][ui] = this->nullWeight;
		}

		// TODO: remove attributes

		m--; // decreasing the number of edges

	}
}

bool Graph::hasEdge(node u, node v) const {
	return (find(u, v) != none);
}


node Graph::mergeEdge(node u, node v, bool discardSelfLoop) {
	DEBUG("merge edge with nodes " , u , " and " , v);

	if (u != v) {
		node newNode = this->addNode();

		// self-loop if necessary
		if (! discardSelfLoop) {
			TRACE("selfLoopWeight");
			edgeweight selfLoopWeight = this->weight(u, u) + this->weight(v, v) + this->weight(u, v);
			this->addEdge(newNode, newNode, selfLoopWeight);
			TRACE("end selfLoopWeight");
		}

		// rewire edges from u to newNode
		this->forWeightedEdgesOf(u, [&](node u, node neighbor, edgeweight w) {
			if (neighbor != u && neighbor != v) {
				TRACE("neighbor of " , u , ": " , neighbor);
				this->increaseWeight(neighbor, newNode, this->weight(u, neighbor)); // TODO: make faster
				TRACE("end neighbor of u");
			}
		});

		// rewire edges from v to newNode
		this->forWeightedEdgesOf(v, [&](node v, node neighbor, edgeweight w) {
			if (neighbor != v && neighbor != u) {
				TRACE("neighbor of " , v , ": " , neighbor);
				this->increaseWeight(neighbor, newNode, this->weight(v, neighbor));  // TODO: make faster
				TRACE("end neighbor of v");
			}
		});

		// delete edges of nodes to delete
		this->forEdgesOf(u, [&](node u, node neighbor) {
			this->removeEdge(u, neighbor);
		});
		this->forEdgesOf(v, [&](node v, node neighbor) {
			this->removeEdge(v, neighbor);
		});

		TRACE("incident edges deleted");

		// delete nodes
		this->removeNode(u);
		this->removeNode(v);

		TRACE("u and v deleted");

		return newNode;
	}

	// no new node created
	return none;
}


/** EDGE ATTRIBUTES **/

edgeweight Graph::weight(node u, node v) const {
	index vi = find(u, v);
	if (vi != none) {
		if (weighted) {
			return this->eweights[u][vi];
		} else {
			return defaultEdgeWeight;
		}
	} else {
		return 0.0;
	}
}

void Graph::setWeight(node u, node v, edgeweight w) {
	if (!weighted) {
		throw std::runtime_error("this is an unweighted graph");	
	}
	if (u == v) {
		// self-loop case
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

void Graph::increaseWeight(node u, node v, edgeweight w) {
	if (!weighted) {
		throw std::runtime_error("this is an unweighted graph");
	}
	if (this->hasEdge(u, v)) {
		this->setWeight(u, v, w + this->weight(u, v));
	} else {
		this->addEdge(u, v, w);
	}
}

int Graph::addEdgeAttribute_double(double defaultValue) {
	int attrId = this->edgeMaps_double.size();

	std::vector<std::vector<double> > edgeMap(z);
	if (this->numberOfEdges() > 0) {
		forNodes([&] (node v) {
			// create edgeMaps and fill them with default value
			edgeMap[v].resize(adja[v].size());
			fill(edgeMap[v].begin(), edgeMap[v].end(), defaultValue);
		});
	}
	
	this->edgeMaps_double.push_back(edgeMap);
	this->edgeAttrDefaults_double.push_back(defaultValue);

	return attrId;
}

double Graph::attribute_double(node u, node v, int attrId) const {
	assert (attrId < (int) this->edgeMaps_double.size());
	index vi = find(u, v);
	if (vi != none) {
		return this->edgeMaps_double[attrId][u][vi];
	} else {
		throw std::runtime_error("What if edge does not exist?");
	}
}

void Graph::setAttribute_double(node u, node v, int attrId, double attr) {
	if (u == v) {
		// self-loop case
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
			this->edgeMaps_double[attrId][u][vi] = attr;
			this->edgeMaps_double[attrId][v][ui] = attr;
		} else {
			throw std::runtime_error("What if edge does not exist?");
		}
	}
}


/** GLOBAL PROPERTIES **/

bool Graph::isDirected() const {
	return false;
}

} /* namespace NetworKit */
