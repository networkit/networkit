/*
 * BasicGraph.cpp
 *
 *  Created on: 20.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#include <stdexcept>
#include <algorithm>
#include <sstream>

#include "BasicGraph.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

namespace graph_impl {

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

template<Weighted w>
index find_impl(const BasicGraph<w, Directed::undirected>& G, node v) {

for (index i = 0; i < G.adja[v].size(); i++) {
		node x = G.adja[v][i];
		if (x == v) {
			return i;
		}
	}
	return none;
}

template<Weighted w>
index find_impl(const BasicGraph<w, Directed::directed>& G, node v) {

	for (index i = 0; i < G.OutEdges[v].size(); i++) {
		node x = G.OutEdges[v];
		if (x == v) {
			return i;
		}
	}

	return none;


}

template<Weighted w>
count degree_impl(const BasicGraph<w, Directed::undirected>& G, node v) {
	return G.deg[v];
}

template<Weighted w>
count degree_impl(const BasicGraph<w, Directed::directed>& G, node v) {
	return G.outDeg[v];
}

/*EDGE-MODIFIERS*/

template<Weighted w>
void add_Edge_impl(const BasicGraph<w, Directed::undirected>& G, node u, node v, edgeweight weight) {
	
	assert (u >= 0);
	assert (u < this->z);
	assert (this->exists[u]);
	assert (v >= 0);
	assert (v < this->z);
	assert (this->exists[v]);

	if (u == v) { // self-loop case
		G.adja[u].push_back(u);
		G.deg[u] += 1;
		if (w == Weighted::weighted) G.weights[u].push_back(weight);
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
			G.weights[u].push_back(weight);
			G.weights[v].push_back(weight);	
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
void add_Edge_impl(const BasicGraph<w, Directed::directed>& G, node u, node v, edgeweight weight) {

assert (u >= 0);
	assert (u < G.z);
	assert (v >= 0);
	assert (v < G.z);

	G.m++; // increasing the number of edges
	G.deg[u].out++;
	G.deg[v].in++;

	G.adja[u].first.push_back(v);
	if (w == Weighted::weighted) {
		G.weights[u].push_back(weight);
	}
	G.adja[v].second.push_back(u);



	// loop over all attributes, setting default attr
	for (index attrId = 0; attrId < G.edgeMaps_double.size(); ++attrId) {
		double defaultAttr = G.edgeAttrDefaults_double[attrId];
		G.edgeMaps_double[attrId][u].push_back(defaultAttr);
		G.edgeMaps_double[attrId][v].push_back(defaultAttr);
	}

}

template<Weighted w>
void remove_Edge_impl(const BasicGraph<w, Directed::undirected>& G, node u, node v, edgeweight weight) {

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
void remove_Edge_impl(const BasicGraph<w, Directed::directed>& G, node u, node v, edgeweight weight) {

}
template<> int BasicGraph<Weighted::unweighted, Directed::undirected>::typ() const { return 1; }
template<> int BasicGraph<Weighted::weighted,   Directed::undirected>::typ() const { return 2; }
template<> int BasicGraph<Weighted::unweighted, Directed::directed>::typ()   const { return 3; }
template<> int BasicGraph<Weighted::weighted,   Directed::directed>::typ()   const { return 4; }

} // namespace graph_impl



// algorithm for all graph classes
template<Weighted w, Directed d>
count degreeSum(const BasicGraph<w, d>& G) {
	count sum = 0;
	G.forNodes([&](node v) {
		sum += G.degree(v);
	});
	return sum;
}

// algorithm for all undirected graph classes
template<Weighted w>
count undirected_algo(const IUndirectedGraph<w>& G) { 
	return 3141592;
}

template<Weighted w>
count directed_algo(const IDirectedGraph<w>& G) {
	return 2718281828;
}




} /* namespace NetworKit */
