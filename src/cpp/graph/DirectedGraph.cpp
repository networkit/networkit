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
	adja(z) {
	if (weighted) eweights.resize(z); // edge weight array is only initialized for weighted graphs
}

count DirectedGraph::getMemoryUsage() const {
	count mem = AbstractGraph::getMemoryUsage();
	mem += sizeof(NodeDegree) * deg.capacity();
		
	mem += sizeof(std::vector<node>) * adja.capacity();
	for (auto& a : adja) {
		mem += sizeof(node) * a.first.capacity();
	}

	for(auto& a: adja) {
		mem+= sizeof(node) * a.second.capacity();
	}
	//mem += sizeof(index) * inOut.capacity();

	mem += sizeof(std::vector<edgeweight>) * eweights.capacity();
	for (auto& w : eweights) {
		mem += sizeof(edgeweight) * w.capacity();
	}

	for (auto& edgeMap : edgeMaps_double) {
		for (auto& em : edgeMap) {
			mem += sizeof(double) * em.capacity();
		}
	}

	return mem;
}

void DirectedGraph::shrinkToFit() {
	AbstractGraph::shrinkToFit();

	deg.shrink_to_fit();
	//inOut.shrink_to_fit();

	adja.shrink_to_fit();
	for (auto& a : adja) {
		a.first.shrink_to_fit();
		a.second.shrink_to_fit();
	}
	

	eweights.shrink_to_fit();
	for (auto& w : eweights) {
		w.shrink_to_fit();
	}

	for (auto& edgeMap : edgeMaps_double) {
		edgeMap.shrink_to_fit();
		for (auto& em : edgeMap) {
			em.shrink_to_fit();
		}
	}
}

//only to be used by Cython
void DirectedGraph::stealFrom(DirectedGraph& input) {
	*this = std::move(input);
}


index DirectedGraph::findIn(node u, node v) const {
	for (index i = 0; i < (adja[u].second).size(); i++) {
		node x = this->adja[u].second[i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

index DirectedGraph::findOut(node u, node v) const {
	for (index i = 0; i < (adja[u].first).size(); i++) {
		node x = this->adja[u].first[i];
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
	this->adja.push_back(std::pair<std::vector<node>,std::vector<node> >{});
	if (weighted) {
		// vector of edge weights for new node
		this->eweights.push_back(std::vector<edgeweight>{});
	}

	// update edge attribute data structures
	for (int attrId = 0; attrId < (int) this->edgeMaps_double.size(); ++attrId) {
		std::vector<double> attrVector;
		this->edgeMaps_double[attrId].push_back(attrVector);
	}

	return v;
}

node DirectedGraph::addNode(float x, float y) {
	node v = addNode();
	std::vector<float> coords = {x, y};
	coordinates.addCoordinates(coords);
	return v;
}


/** NODE PROPERTIES **/

count DirectedGraph::degreeIn(node v) const {
	assert (v >= 0);
	assert (v < this->z);
	assert (this->exists[v]);
	return this->deg[v].in;
}

count DirectedGraph::degreeOut(node v) const {
	assert (v >= 0);
	assert (v < this->z);
	assert (this->exists[v]);
	return this->deg[v].out;
}

edgeweight DirectedGraph::weightedDegreeIn(node v) const {
	if (!isWeighted()) {
		return degreeOut(v);
	}
	edgeweight sum = 0.0;
	forWeightedOutNeighborsOf(v, [&](node u, edgeweight ew) {
		sum += ew;
	});
	return sum;
}

edgeweight DirectedGraph::weightedDegreeOut(node v) const {
	if (!isWeighted()) {
		return degreeIn(v);
	}
	edgeweight sum = 0.0;
	forWeightedInNeighborsOf(v, [&](node u, edgeweight ew) {
		sum += ew;
	});
	return sum;
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

	this->adja[u].first.push_back(v);
	if (weighted) {
		this->eweights[u].push_back(weight);
	}
	this->adja[v].second.push_back(u);



	// loop over all attributes, setting default attr
	for (index attrId = 0; attrId < this->edgeMaps_double.size(); ++attrId) {
		double defaultAttr = this->edgeAttrDefaults_double[attrId];
		this->edgeMaps_double[attrId][u].push_back(defaultAttr);
		this->edgeMaps_double[attrId][v].push_back(defaultAttr);
	}
}

void DirectedGraph::removeEdge(node u, node v) {
	// remove adjacency
	index ui = findIn(v, u);
	index vi = findOut(u, v);

	if (vi == none) {
		std::stringstream strm;
		strm << "edge (" << u << "," << v << ") does not exist";
		throw std::runtime_error(strm.str());
	} else {
		this->adja[u].first[vi] = none;
		this->adja[v].second[ui] = none;
		// decrement degree counters
		this->deg[u].out--;
		this->deg[v].in--;
		if (weighted) {
			// remove edge weight
			this->eweights[u][vi] = this->nullWeight;
			this->eweights[v][ui] = this->nullWeight;
		}

		// TODO: remove attributes

		m--; // decreasing the number of edges

	}
}

bool DirectedGraph::hasEdge(node u, node v) const {
	return (findOut(u, v) != none);
}



/** EDGE ATTRIBUTES **/

edgeweight DirectedGraph::weight(node u, node v) const {
	index vi = findOut(u, v);
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

void DirectedGraph::setWeight(node u, node v, edgeweight w) {
	if (!weighted) {
		throw std::runtime_error("this is an unweighted graph");	
	}

	index vi = findOut(u, v);
	index ui = findIn(v, u);
	if ((vi != none) && (ui != none)) {
		this->eweights[u][vi] = w;
		this->eweights[v][ui] = w;
	} else {
		addEdge(u, v, w);
	}
}

void DirectedGraph::increaseWeight(node u, node v, edgeweight w) {
	if (!weighted) {
		throw std::runtime_error("this is an unweighted graph");
	}
	if (this->hasEdge(u, v)) {
		this->setWeight(u, v, w + this->weight(u, v));
	} else {
		this->addEdge(u, v, w);
	}
}

int DirectedGraph::addEdgeAttribute_double(double defaultValue) {
	int attrId = this->edgeMaps_double.size();

	std::vector<std::vector<double> > edgeMap(z);
	if (this->numberOfEdges() > 0) {
		forNodes([&] (node v) {
			// create edgeMaps and fill them with default value
			edgeMap[v].resize(adja[v].first.size());
			fill(edgeMap[v].begin(), edgeMap[v].end(), defaultValue);
		});
	}
	
	this->edgeMaps_double.push_back(edgeMap);
	this->edgeAttrDefaults_double.push_back(defaultValue);

	return attrId;
}

double DirectedGraph::attribute_double(node u, node v, int attrId) const {
	assert (attrId < (int) this->edgeMaps_double.size());
	index vi = findOut(u, v);
	if (vi != none) {
		return this->edgeMaps_double[attrId][u][vi];
	} else {
		throw std::runtime_error("What if edge does not exist?");
	}
}

void DirectedGraph::setAttribute_double(node u, node v, int attrId, double attr) {
	index vi = findOut(u, v);
	index ui = findIn(v, u);
	if ((vi != none) && (ui != none)) {
		this->edgeMaps_double[attrId][u][vi] = attr;
		this->edgeMaps_double[attrId][v][ui] = attr;
	} else {
		throw std::runtime_error("What if edge does not exist?");
	}
}


/** GLOBAL PROPERTIES **/

bool DirectedGraph::isDirected() const {
	return true;
}


} /* namespace NetworKit */
