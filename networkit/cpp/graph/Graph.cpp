/*
 * Graph.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#include <sstream>
#include <random>

#include "Graph.h"
#include "GraphBuilder.h"

namespace NetworKit {

/** CONSTRUCTORS **/

Graph::Graph(count n, bool weighted, bool directed) :
	n(n),
	m(0),
	storedNumberOfSelfLoops(0),
	z(n),
	omega(0),
	t(0),

	weighted(weighted), // indicates whether the graph is weighted or not
	directed(directed), // indicates whether the graph is directed or not
	edgesIndexed(false), // edges are not indexed by default

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
	outEdgeWeights(weighted ? n : 0),
	inEdgeIds(),
	outEdgeIds() {

	// set name from global id
	id = getNextGraphId();
	std::stringstream sstm;
	sstm << "G#" << id;
	name = sstm.str();
}

Graph::Graph(std::initializer_list<WeightedEdge> edges) : Graph(0, true) {
  using namespace std;

  /* Number of nodes = highest node index + 1 */
  for (const auto& edge: edges) {
    node x = max(edge.u, edge.v);
    while (numberOfNodes() <= x) {
      addNode();
    }
  }

  /* Now add all of the edges */
  for (const auto& edge: edges) {
    addEdge(edge.u, edge.v, edge.weight);
  }
}

Graph::Graph(const Graph& G, bool weighted, bool directed) :
	n(G.n),
	m(G.m),
	storedNumberOfSelfLoops(G.storedNumberOfSelfLoops),
	z(G.z),
	omega(0),
	t(G.t),
	weighted(weighted),
	directed(directed),
	edgesIndexed(false), //edges are not indexed by default
	exists(G.exists),

	// let the following be empty for the start, we fill them later
	inDeg(0),
	outDeg(0),
	inEdges(0),
	outEdges(0),
	inEdgeWeights(0),
	outEdgeWeights(0) {

	// set name from global id
	id = getNextGraphId();
	std::stringstream sstm;
	sstm << "G#" << id;
	name = sstm.str();

	if (G.isDirected() == directed) {
		inDeg = G.inDeg; // G.inDeg might be empty (if G is undirected), but that's fine
		outDeg = G.outDeg;
		inEdges = G.inEdges; // G.inEdges might be empty (if G is undirected), but that's fine
		outEdges = G.outEdges;

		// copy weights if needed
		if (weighted) {
			if (G.isWeighted()) {
				// just copy from G, again either both graphs are directed or both are undirected
				inEdgeWeights = G.inEdgeWeights;
				outEdgeWeights = G.outEdgeWeights;
			} else {
				// G has no weights, set defaultEdgeWeight for all edges
				if (directed) {
					inEdgeWeights.resize(z);
					for (node u = 0; u < z; u++) {
						inEdgeWeights[u] = std::vector<edgeweight>(G.inEdges[u].size(), defaultEdgeWeight);
					}
				}

				outEdgeWeights.resize(z);
				for (node u = 0; u < z; u++) {
					outEdgeWeights[u] = std::vector<edgeweight>(outEdges[u].size(), defaultEdgeWeight);
				}
			}
		}
	} else if (G.isDirected()) {
		// G is directed, but we want an undirected graph
		// so we need to combine the out and in stuff for every node
		outDeg.resize(z);
		outEdges.resize(z);
		for (node u = 0; u < z; u++) {
			outDeg[u] = G.inDeg[u] + G.outDeg[u];

			// copy both out and in edges into our new outEdges
			outEdges[u].reserve(G.outEdges[u].size() + G.inEdges[u].size());
			outEdges[u].insert(outEdges[u].end(), G.outEdges[u].begin(), G.outEdges[u].end());
			outEdges[u].insert(outEdges[u].end(), G.inEdges[u].begin(), G.inEdges[u].end());
		}
		if (weighted) {
			if (G.isWeighted()) {
				// same for weights
				outEdgeWeights.resize(z);
				for (node u = 0; u < z; u++) {
					outEdgeWeights[u].reserve(G.outEdgeWeights[u].size() + G.inEdgeWeights[u].size());
					outEdgeWeights[u].insert(outEdgeWeights[u].end(), G.outEdgeWeights[u].begin(), G.outEdgeWeights[u].end());
					outEdgeWeights[u].insert(outEdgeWeights[u].end(), G.inEdgeWeights[u].begin(), G.inEdgeWeights[u].end());
				}
			} else {
				// we are undirected, so no need to write anything into inEdgeWeights
				outEdgeWeights.resize(z);
				for (node u = 0; u < z; u++) {
					outEdgeWeights[u] = std::vector<edgeweight>(outEdges[u].size(), defaultEdgeWeight);
				}
			}
		}
	} else {
		// G is not directed, but this copy should be
		// generally we can can copy G.out stuff into our in stuff
		inDeg = G.outDeg;
		outDeg = G.outDeg;
		inEdges = G.outEdges;
		outEdges = G.outEdges;
		if (weighted) {
			if (G.isWeighted()) {
				inEdgeWeights = G.outEdgeWeights;
				outEdgeWeights = G.outEdgeWeights;
			} else {
				// initialize both inEdgeWeights and outEdgeWeights with the defaultEdgeWeight
				inEdgeWeights.resize(z);
				for (node u = 0; u < z; u++) {
					inEdgeWeights[u] = std::vector<edgeweight>(inEdges[u].size(), defaultEdgeWeight);
				}
				outEdgeWeights.resize(z);
				for (node u = 0; u < z; u++) {
					outEdgeWeights[u] = std::vector<edgeweight>(outEdges[u].size(), defaultEdgeWeight);
				}
			}
		}
	}

}

/** PRIVATE HELPERS **/

count Graph::getNextGraphId() {
	static count nextGraphId = 1;
	return nextGraphId++;
}

index Graph::indexInInEdgeArray(node v, node u) const {
	if (!directed) {
		return indexInOutEdgeArray(v, u);
	}
	for (index i = 0; i < inEdges[v].size(); i++) {
		node x = inEdges[v][i];
		if (x == u) {
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


/** EDGE IDS **/

void Graph::indexEdges(bool force) {
	if (edgesIndexed && !force) return;

	omega = 0; // reset edge ids (for re-indexing)

	outEdgeIds.resize(outEdges.size());
	forNodes([&](node u) {
		outEdgeIds[u].resize(outEdges[u].size(), none);
	});

	if (directed) {
		inEdgeIds.resize(inEdges.size());
		forNodes([&](node u) {
			inEdgeIds[u].resize(inEdges[u].size(), none);
		});
	}

	// assign edge ids for edges in one direction
	forNodes([&](node u) {
		for (index i = 0; i < outEdges[u].size(); ++i) {
			node v = outEdges[u][i];
			if (v != none && (directed || (u >= v))) {
				// new id
				edgeid id = omega++;
				outEdgeIds[u][i] = id;
			}
		}
	});

	// copy edge ids for the edges in the other direction. Note that "indexInOutEdgeArray" is slow
	// which is why this second loop in parallel makes sense.
	if (!directed) {
		balancedParallelForNodes([&](node u) {
			for (index i = 0; i < outEdges[u].size(); ++i) {
				node v = outEdges[u][i];
				if (v != none && outEdgeIds[u][i] == none) {
					index j = indexInOutEdgeArray(v, u);
					outEdgeIds[u][i] = outEdgeIds[v][j];
				}
			}
		});
	} else {
		balancedParallelForNodes([&](node u) {
			for (index i = 0; i < inEdges[u].size(); ++i) {
				node v = inEdges[u][i];
				if (v != none) {
					index j = indexInOutEdgeArray(v, u);
					inEdgeIds[u][i] = outEdgeIds[v][j];
				}
			}
		});
	}

	edgesIndexed = true; // remember that edges have been indexed so that addEdge needs to create edge ids
}


edgeid Graph::edgeId(node u, node v) const {
	if (! edgesIndexed) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}
	index i = indexInOutEdgeArray(u, v);

	if (i == none) {
		throw std::runtime_error("Edge does not exist");
	}

	edgeid id = outEdgeIds[u][i];
	return id;
}


/** GRAPH INFORMATION **/

std::string Graph::typ() const {
	if (weighted) {
		return directed ? "WeightedDirectedGraph" : "WeightedGraph";
	} else {
		return directed ? "DirectedGraph" : "Graph";
	}
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

}

void Graph::compactEdges() {
	this->parallelForNodes([&](node u) {
		if (degreeOut(u) != outEdges[u].size()) {
			if (degreeOut(u) == 0) {
				outEdges[u].clear();
				if (weighted) outEdgeWeights[u].clear();
				if (edgesIndexed) outEdgeIds[u].clear();
			} else {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					while (i < outEdges[u].size() && outEdges[u][i] == none) {
						outEdges[u][i] = outEdges[u].back();
						outEdges[u].pop_back();

						if (weighted) {
							outEdgeWeights[u][i] = outEdgeWeights[u].back();
							outEdgeWeights[u].pop_back();
						}

						if (edgesIndexed) {
							outEdgeIds[u][i] = outEdgeIds[u].back();
							outEdgeIds[u].pop_back();
						}
					}
				}
			}
		}

		if (directed && degreeIn(u) != inEdges[u].size()) {
			if (degreeIn(u) == 0) {
				inEdges[u].clear();
				if (weighted) inEdgeWeights[u].clear();
			       if (edgesIndexed) inEdgeIds[u].clear();
			} else {
				for (index i = 0; i < inEdges[u].size(); ++i) {
					while (i < inEdges[u].size() && inEdges[u][i] == none) {
						inEdges[u][i] = inEdges[u].back();
						inEdges[u].pop_back();

						if (weighted) {
							inEdgeWeights[u][i] = inEdgeWeights[u].back();
							inEdgeWeights[u].pop_back();
						}

						if (edgesIndexed) {
							inEdgeIds[u][i] = inEdgeIds[u].back();
							inEdgeIds[u].pop_back();
						}
					}
				}
			}

		}
	});
}

void Graph::sortEdges() {
	std::vector<std::vector<node> > targetAdjacencies(upperNodeIdBound());
	std::vector<std::vector<edgeweight> > targetWeight;
	std::vector<std::vector<edgeid> > targetEdgeIds;

	if (isWeighted()) {
		targetWeight.resize(upperNodeIdBound());
		forNodes([&](node u) {
			targetWeight[u].reserve(degree(u));
		});
	}
	if (hasEdgeIds()) {
		targetEdgeIds.resize(upperNodeIdBound());
		forNodes([&](node u) {
			targetEdgeIds[u].reserve(degree(u));
		});
	}

	forNodes([&](node u) {
		targetAdjacencies[u].reserve(degree(u));
	});

	auto assignToTarget = [&](node u, node v, edgeweight w, edgeid eid) {
			targetAdjacencies[v].push_back(u);
			if (isWeighted()) {
				targetWeight[v].push_back(w);
			}
			if (hasEdgeIds()) {
				targetEdgeIds[v].push_back(eid);
			}
		};

	forNodes([&](node u) {
		forInEdgesOf(u, assignToTarget);
	});

	outEdges.swap(targetAdjacencies);
	outEdgeWeights.swap(targetWeight);
	outEdgeIds.swap(targetEdgeIds);

	if (isDirected()) {
		inEdges.swap(targetAdjacencies);
		inEdgeWeights.swap(targetWeight);
		inEdgeIds.swap(targetEdgeIds);

		forNodes([&](node u) {
			targetAdjacencies[u].resize(degreeIn(u));
			targetAdjacencies[u].shrink_to_fit();
			targetAdjacencies[u].clear();
			if (isWeighted()) {
				targetWeight[u].resize(degreeIn(u));
				targetWeight[u].shrink_to_fit();
				targetWeight[u].clear();
			}
			if (hasEdgeIds()) {
				targetEdgeIds[u].resize(degreeIn(u));
				targetEdgeIds[u].shrink_to_fit();
				targetEdgeIds[u].clear();
			}
		});

		forNodes([&](node u) {
			forEdgesOf(u, assignToTarget);
		});

		inEdges.swap(targetAdjacencies);
		inEdgeWeights.swap(targetWeight);
		inEdgeIds.swap(targetEdgeIds);
	}
}


std::string Graph::toString() const {
	std::stringstream strm;
	strm << typ() << "(name=" << getName() << ", n=" << numberOfNodes()
			<< ", m=" << numberOfEdges() << ")";
	return strm.str();
}


/** COPYING **/

Graph Graph::copyNodes() const {
	Graph C(z, weighted, directed);
	for (node u = 0; u < z; ++u) {
		if (! exists[u]) {
			C.removeNode(u);
		}
	}
	return C;
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

void Graph::restoreNode(node v){
	assert(v < z);
	assert(!exists[v]);

	exists[v] = true;
	n++;
}


/** NODE PROPERTIES **/

edgeweight Graph::weightedDegree(node v) const {
	if (weighted) {
		edgeweight sum = 0.0;
		forNeighborsOf(v, [&](node u, edgeweight ew) {
			sum += ew;
		});
		return sum;
	}
	return defaultEdgeWeight * degree(v);
}

edgeweight Graph::volume(node v) const {
	if (weighted) {
		edgeweight sum = 0.0;
		for (index i = 0; i < outEdges[v].size(); i++) {
			node u = outEdges[v][i];
			if (u == v) {
				sum += 2 * outEdgeWeights[v][i];
			} else if (u != none) {
				sum += outEdgeWeights[v][i];
			}
		}
		return sum;
	} else {
		count c = outDeg[v];
		for (node u : outEdges[v]) {
			if (u == v) {
				c++;
			}
		}
		return c * defaultEdgeWeight;
	}
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
	if (outDeg[u] == 0) {
		return none;
	}

	node v;
	do {
		index i = Aux::Random::integer(outEdges[u].size() - 1);
		v = outEdges[u][i];
	} while (v == none);

	return v;
}


/** EDGE MODIFIERS **/

void Graph::addEdge(node u, node v, edgeweight ew) {
	assert (u < z);
	assert (exists[u]);
	assert (v < z);
	assert (exists[v]);

	m++; // increase number of edges
	outDeg[u]++;
	outEdges[u].push_back(v);

	// if edges indexed, give new id
	if (edgesIndexed) {
		edgeid id = omega++;
		outEdgeIds[u].push_back(id);
	}

	if (directed) {
		inDeg[v]++;
		inEdges[v].push_back(u);

		if (edgesIndexed) {
			inEdgeIds[v].push_back(id);
		}

		if (weighted) {
			inEdgeWeights[v].push_back(ew);
			outEdgeWeights[u].push_back(ew);
		}

	} else if (u == v) { // self-loop case
		if (weighted) {
			outEdgeWeights[u].push_back(ew);
		}

	} else { // undirected, no self-loop
		outDeg[v]++;
		outEdges[v].push_back(u);

		if (weighted) {
			outEdgeWeights[u].push_back(ew);
			outEdgeWeights[v].push_back(ew);
		}

		if (edgesIndexed) {
			outEdgeIds[v].push_back(omega - 1);
		}
	}

	if (u == v) { //count self loop
		storedNumberOfSelfLoops++;
	}
}

void Graph::removeEdge(node u, node v) {
	assert (u < z);
	assert (exists[u]);
	assert (v < z);
	assert (exists[v]);

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
		assert (ui != none);

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

	if (u == v) {
		storedNumberOfSelfLoops--;
		assert(storedNumberOfSelfLoops >= 0);
	}

	// dose not make a lot of sense do remove attributes,
	// cause the edge is marked as deleted and we have no null values for the attributes
}

void Graph::removeSelfLoops() {
	this->forEdges([&](node u, node v, edgeweight ew) {
		if (u == v) {
			removeEdge(u, v);
		}
	});

}


void Graph::swapEdge(node s1, node t1, node s2, node t2) {
	index s1t1 = indexInOutEdgeArray(s1, t1);
	if (s1t1 == none) throw std::runtime_error("The first edge does not exist");
	index t1s1 = indexInInEdgeArray(t1, s1);

	index s2t2 = indexInOutEdgeArray(s2, t2);
	if (s2t2 == none) throw std::runtime_error("The second edge does not exist");
	index t2s2 = indexInInEdgeArray(t2, s2);

	std::swap(outEdges[s1][s1t1], outEdges[s2][s2t2]);

	if (directed) {
		std::swap(inEdges[t1][t1s1], inEdges[t2][t2s2]);

		if (weighted) {
			std::swap(inEdgeWeights[t1][t1s1], inEdgeWeights[t2][t2s2]);
		}

		if (edgesIndexed) {
			std::swap(inEdgeIds[t1][t1s1], inEdgeIds[t2][t2s2]);
		}
	} else {
		std::swap(outEdges[t1][t1s1], outEdges[t2][t2s2]);

		if (weighted) {
			std::swap(outEdgeWeights[t1][t1s1], outEdgeWeights[t2][t2s2]);
		}

		if (edgesIndexed) {
			std::swap(outEdgeIds[t1][t1s1], outEdgeIds[t2][t2s2]);
		}
	}
}

bool Graph::hasEdge(node u, node v) const {
	return indexInOutEdgeArray(u, v) != none;
}

std::pair<node, node> Graph::randomEdge(bool uniformDistribution) const {
	if (m == 0) {
		return std::make_pair(none, none);
	}

	if (uniformDistribution) {
		return randomEdges(1)[0];
	}

	node u, v; // we will return edge (u, v)
	// fast way, but not a uniform random edge!
	do {
		u = randomNode();
	} while (outDeg[u] == 0);
	v = randomNeighbor(u);
	return std::make_pair(u, v);
}

std::vector< std::pair<node, node> > Graph::randomEdges(count nr) const {
	std::vector< std::pair<node, node> > edges;

	std::default_random_engine gen{std::random_device{}()};
	std::discrete_distribution<count> distribution(outDeg.begin(), outDeg.end());

	for (index i = 0; i < nr; i++) {
		node u, v; // we will pick edge (u, v)
		if (directed) {
			u = distribution(gen);
			assert(outEdges[u].size() > 0); // should always be the case as  without edges should have probability 0
			v = randomNeighbor(u);
		} else {
			// self-loops which appear only once in the outEdge arrays
			// easiest way it to ignore edges (u, v) with u > v
			do {
				u = distribution(gen);
				assert(outEdges[u].size() > 0); // should always be the case as  without edges should have probability 0
				v = randomNeighbor(u);
			} while (u > v);
		}
		edges.push_back({u, v});
	}

	return edges;
}


/** GLOBAL PROPERTIES **/

count Graph::numberOfSelfLoops() const {
	return storedNumberOfSelfLoops;
}


/** EDGE ATTRIBUTES **/

edgeweight Graph::weight(node u, node v) const {
	index vi = indexInOutEdgeArray(u, v);
	if (vi == none) {
		return nullWeight;
	} else {
		return weighted ? outEdgeWeights[u][vi] : defaultEdgeWeight;
	}
}

void Graph::setWeight(node u, node v, edgeweight ew) {
	if (!weighted) {
		throw std::runtime_error("Cannot set edge weight in unweighted graph.");
	}

	index vi = indexInOutEdgeArray(u, v);
	if (vi == none) {
		// edge does not exist, create it, but warn user
		TRACE("Setting edge weight of a nonexisting edge will create the edge.");
		addEdge(u, v, ew);
		return;
	}

	outEdgeWeights[u][vi] = ew;
	if (directed) {
		index ui = indexInInEdgeArray(v, u);
		inEdgeWeights[v][ui] = ew;
	} else if (u != v) {
		index ui = indexInInEdgeArray(v, u);
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
		TRACE("Increasing edge weight of a nonexisting edge will create the edge.");
		addEdge(u, v, ew);
		return;
	}

	outEdgeWeights[u][vi] += ew;
	if (directed) {
		index ui = indexInInEdgeArray(v, u);
		inEdgeWeights[v][ui] += ew;
	} else if (u != v) {
		index ui = indexInInEdgeArray(v, u);
		outEdgeWeights[v][ui] += ew;
	}
}



/** SUMS **/

edgeweight Graph::totalEdgeWeight() const {
	if (weighted) {
		edgeweight sum = 0.0;
		forEdges([&](node u, node v, edgeweight ew) {
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

Graph Graph::transpose() const {
	if (directed == false) {
		throw std::runtime_error("The transpose of an undirected graph is identical to the original graph.");
	}

	GraphBuilder gB(z, weighted, true);

	this->forEdges([&](node u, node v) {
		gB.addHalfEdge(v, u, weight(u,v));
	});

	Graph GTranspose = gB.toGraph(true);
	for (node u = 0; u < z; ++u) {
		if (! exists[u]) {
			GTranspose.removeNode(u);
		}
	}
	GTranspose.t = t;
	GTranspose.setName(getName() + "Transpose");
	return GTranspose;
}

Graph Graph::toUndirected() const {
	if (directed == false) {
		throw std::runtime_error("this graph is already undirected");
	}
	Graph U(*this, weighted, false);
	return U;
}


Graph Graph::toUnweighted() const {
	if (weighted == false) {
		throw std::runtime_error("this graph is already unweighted");
	}
	Graph U(*this, false, directed);
	return U;
}

bool Graph::checkConsistency() const {
	// check for multi-edges
	std::vector<node> lastSeen(z, none);
	bool noMultiEdges = true;
	auto noMultiEdgesDetected = [&noMultiEdges]() { return noMultiEdges; };
	forNodesWhile(noMultiEdgesDetected, [&](node v) {
		forNeighborsOf(v, [&](node u) {
			if (lastSeen[u] == v) {
				noMultiEdges = false;
				DEBUG("Multiedge found between ", u, " and ", v, "!");
			}
			lastSeen[u] = v;
		});
	});

	return noMultiEdges;
}

void Graph::append(const Graph& G) {
	std::map<node,node> nodeMap;
	G.forNodes([&](node u) {
		node u_ = this->addNode();
		nodeMap[u] = u_;
	});
	if (this->isWeighted()) {
		G.forEdges([&](node u, node v, edgeweight ew){
			this->addEdge(nodeMap[u], nodeMap[v], ew);
		});
	} else {
		G.forEdges([&](node u, node v){
			this->addEdge(nodeMap[u], nodeMap[v]);
		});
	}
}

void Graph::merge(const Graph& G) {
	// TODO: handle edge weights
	G.forEdges([&](node u, node v) {
		// naive implementation takes $O(m \cdot d)$ for $m$ edges and max. degree $d$ in this graph
		if (!this->hasEdge(u, v)) {
			this->addEdge(u, v);
		}
	});
}


// SUBGRAPHS


Graph Graph::subgraphFromNodes(const std::unordered_set<node>& nodes) const {

	Graph S(upperNodeIdBound(), isWeighted(), isDirected());
	// delete all nodes that are not in the node set
	S.forNodes([&](node u) {
		if (nodes.find(u) == nodes.end()) {
			S.removeNode(u);
		}
	});

	forEdges([&](node u, node v, edgeweight w) {
		// if both end nodes are in the node set
		if (nodes.find(u) != nodes.end() && nodes.find(v) != nodes.end()) {
			S.addEdge(u, v, w);
		}
	});

	return S;
}



} /* namespace NetworKit */
