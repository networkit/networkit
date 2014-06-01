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
count degree_impl(const BasicGraph<w, Directed::undirected>& G, node v) {
	return G.deg[v];
}

template<Weighted w>
count degree_impl(const BasicGraph<w, Directed::directed>& G, node v) {
	return G.outDeg[v];
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