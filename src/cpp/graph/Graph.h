/*
 * BasicGraph.h
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef BASICGRAPH_H_
#define BASICGRAPH_H_

#include <algorithm>
#include <vector>
#include <stack>
#include <queue>
#include <utility>
#include <stdexcept>

#include "../Globals.h"
#include "Coordinates.h"
#include "../viz/Point.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

class Graph final {
private:
	// graph attributes
	count id;
	std::string name;

	// scalars
	count n; //!< current number of nodes
	count m; //!< current number of edges
	node z; //!< current upper bound of node ids
	count t; //!< current time step

	bool weighted;
	bool directed;

	// per node data
	std::vector<bool> exists; //!< exists[v] is true if node v has not been removed from the graph
	Coordinates<float> coordinates; //!< coordinates of nodes (if present)

	std::vector<count> inDeg;
	std::vector<count> outDeg;
	
	std::vector< std::vector<node> > inEdges;
	std::vector< std::vector<node> > outEdges;
	
	std::vector< std::vector<edgeweight> > inEdgeWeights;
	std::vector< std::vector<edgeweight> > outEdgeWeights;

	// user-defined edge attributes
	// attribute maps storage
	std::vector<std::vector<std::vector<double> > > edgeMaps_double; // contains edge maps (u, v) -> double
	// default values
	std::vector<double> edgeAttrDefaults_double; // stores default value for edgeMaps_double[i] at index i
	// TODO directed case?

	index indexInInEdgeArray(node u, node v) const;
	index indexInOutEdgeArray(node u, node v) const;

public:

	Graph(count n = 0, bool weighted = false, bool directed = false);

	Graph(const Graph& other) = default;

	Graph(Graph&& other) = default;

	~Graph() = default;

	Graph& operator=(Graph&& other) = default;

	Graph& operator=(const Graph& other) = default;

	/** Only to be used from Cython */
	void stealFrom(Graph& input);


	/** GRAPH INFORMATION **/

	/**
	 * Get the ID of this graph. The ID is a unique unsigned integer given to
	 * every graph on construction.
	 */
	count getId() const { return id; }

	/**
	 * Return the type of the graph.
	 * 		Graph: not weighted, undirected
	 * 		WeightedGraph: weighted, undirected
	 * 		DirectedGraph: not weighted, directed
	 * 		WeightedDirectedGraph: weighted, directed
	 */
	std::string typ() const;

	/**
	 * Calculate an approximation of the memory used by this graph. Only memory increasing with the
	 * number of edges or nodes of this graph is taken into account. 
	 */
	count getApproximatedMemoryUsage() const;

	/**
	 * Try to save some memory by shrinking internal data structures of the graph. Only run this
	 * once you finished editing the graph. Otherwise it will cause unnecessary reallocation of
	 * memory. 
	 */
	void shrinkToFit();

	/**
	 * Set name of graph.
	 */
	void setName(std::string name) { this->name = name; }

	/*
	 * @return name of graph
	 */
	std::string getName() const { return name; }

	/**
	 * Get string representation
	 */
	std::string toString() const;

	/** NODE MODIFIERS **/

	/**
	 * Add a new node to the graph and return it.
	 */
	node addNode();

	/**
	 * Add a new node to the graph with coordinates @a x and @y and return it.
	 */
	node addNode(float x, float y);

	/**
	 * Remove an isolated node v from the graph.
	 *
	 * Although it would be convenient to remove all incident edges at the same time,
	 * this causes complications for dynamic applications. Therefore, removeNode is an
	 * atomic event. All incident edges need to be removed first and an exception is thrown
	 * otherwise.
	 */
	void removeNode(node v);

	/**
	 * Check if node v exists in the graph.
	 */
	bool hasNode(node v) const { return (v < z) && this->exists[v];	}


	/** NODE PROPERTIES **/

	/**
	 * Return the number of neighbors for node v.
	 */
	count degree(node v) const { return outDeg[v]; }
	count degreeIn(node v) const { return directed ? inDeg[v] : outDeg[v]; }
	count degreeOut(node v) const { return outDeg[v]; }

	/**
	 * @return true if the node is isolated (= degree is 0)
	 */
	bool isIsolated(node v) const { return outDeg[v] == 0 && (!directed || inDeg[v] == 0); }

	/**
	 * @return Weighted degree of @a v. For directed graphs this is the sum of weights off all outgoing edges fo @a v.
	 */
	edgeweight weightedDegree(node v) const;

	/**
	 * @return Volume of the node, which is the
	 * weighted degree with self-loops counted twice.
	 */
	edgeweight volume(node v) const;

	/**
	 * @return random node of the graph
	 */
	node randomNode() const;

	// TODO comment
	node randomNeighbor(node u) const;


	/** EDGE MODIFIERS **/

	/**
	 * Insert an directed edge between from @a u to @a v.
	 */
	void addEdge(node u, node v, edgeweight ew = defaultEdgeWeight);

	/**
	 * Remove directed edge between from @a u to @a v.
	 */
	void removeEdge(node u, node v);

	/**
	 * Check if directed edge {u,v} exists.
	 */
	bool hasEdge(node u, node v) const;

	/**
	 * Merges edge {u,v} to become a supernode. Edges to u and v are
	 * rewired, multiple edges merged and their weights added.
	 * The vertex weights of @a u and @a v are added.
	 * A self-loop is only created if @a discardSelfLoop is set to false.
	 *
	 * @return New node that has been created if u != v. Otherwise none.
	 */

	/**
	 * @return Random edge
	 */
	std::pair<node, node> randomEdge() const;

	/** GLOBAL PROPERTIES **/

	/**
	 * Return true if this graph supports edge weights other than 1.0
	 */
	bool isWeighted() const { return weighted; }

	/** 
	 * Return true if this graph supports directed edges.
	 */
	bool isDirected() const { return directed; }

	/**
	 * Return true if graph contains no nodes.
	 */
	bool isEmpty() const { return n == 0; }

	/**
	 * Return the number of nodes in the graph.
	 *
	 */
	count numberOfNodes() const { return n; }

	/**
	 * Return the number of edges in the graph.
	 */
	count numberOfEdges() const { return m; }

	/**
	 * @return the number of loops {v, v} in the graph.
	 *
	 * This involves calculation, so store result if needed multiple times.
	 */
	count numberOfSelfLoops() const;

 	/**
	 * Get an upper bound for the node ids in the graph.
	 */
	index upperNodeIdBound() const { return z; }

	/** DYNAMICS **/

	/**
	 * Trigger a time step - increments counter.
	 */
	void timeStep() { t++; }

	/**
	 * Get time step counter.
	 */
	count time() { return t; }


	/** COORDINATES **/

	void setCoordinate(node v, Point<float> value) { coordinates.setCoordinate(v, value); } 

	Point<float>& getCoordinate(node v) { return coordinates.getCoordinate(v); } 

	float minCoordinate(count dim) { return coordinates.minCoordinate(dim); }

	float maxCoordinate(count dim) { return coordinates.maxCoordinate(dim); }

	void initCoordinates() { coordinates.init(z); }


	/** EDGE ATTRIBUTES **/

	/**
	 * Return edge weight.
	 *
	 * Return 0 if edge does not exist.
	 */
	edgeweight weight(node u, node v) const;

	/**
	 * Set the weight of an edge. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	void setWeight(node u, node v, edgeweight ew);

	/**
	 * Increase the weight of an edge. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	void increaseWeight(node u, node v, edgeweight ew);

	/**
	 * Add new edge map for an attribute of type double.
	 */
	int addEdgeAttribute_double(double defaultValue);

	/**
	 * @return attribute of type double for an edge.
	 *
	 * @param[in]	u	node
	 * @param[in]	v	node
	 * @param[in]	attrId	attribute id
	 */
	double attribute_double(node u, node v, int attrId) const;

	/**
	 * Set edge attribute of type double If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	attr	double edge attribute
	 */
	void setAttribute_double(node u, node v, int attrId, double attr);


	/** SUMS **/

	/**
	 * @return sum of all edge weights
	 */
	edgeweight totalEdgeWeight() const;


	/** Collections **/

	/**
	 * Return list of nodes
	 */
	std::vector<node> nodes() const;

	/**
	 * Return list of edges as node pairs.
	 */
	std::vector<std::pair<node, node> > edges() const;

	/**
	 * Return list of neighbors for given node.
	 */
	std::vector<node> neighbors(node u) const;


	/** NODE ITERATORS **/

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodes(L handle) const;

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodes(L handle) const;

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure) as long as the condition remains true.
	 * This allows for breaking from a node loop.
	 */
	template<typename C, typename L> void forNodesWhile(C condition, L handle) const;

	/**
	 * Iterate randomly over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodesInRandomOrder(L handle) const;

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 * Using schedule(guided) to remedy load-imbalances due to e.g. unequal degree distribution.
	 */
	template<typename L> void balancedParallelForNodes(L handle) const;

	/**
	 * Iterate over all undirected pairs of nodesand call handler (lambda closure).
	 */
	template<typename L> void forNodePairs(L handle) const;

	/**
	 * Iterate over all undirected pairs of nodes in parallel and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodePairs(L handle) const;


	/** EDGE ITERATORS **/

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void forEdges(L handle) const;

	/**
	 * Iterate in parallel over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForEdges(L handle) const;

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 *
	 * Handler takes arguments (u, v, w) where u and v are the nodes of the edge and w is its weight.
	 */
	template<typename L> void forWeightedEdges(L handle) const;

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 *
	 * Handler takes arguments (u, v, w) where u and v are the nodes of the edge and w is its weight.
	 */
	template<typename L> void parallelForWeightedEdges(L handle) const;

	/**
	 * Iterate over all edges of the const graph and call handler (lambda closure).
	 *
	 *	@param[in]	attrId		attribute id
	 *	@param[in]	handle 		takes arguments (u, v, a) where a is an edge attribute of edge {u, v}
	 *
	 */
	template<typename L> void forEdgesWithAttribute_double(int attrId, L handle) const;

	/** NEIGHBORHOOD ITERATORS **/

	/**
	 * Iterate over all neighbors of a node and call handler (lamdba closure).
	 * For directed graphs only outgoing edges from u are considered.
	 */
	template<typename L> void forNeighborsOf(node u, L handle) const;

	/**
	 * Iterate over all edge weights of a node and call handler (lamdba closure).
	 * For directed graphs only outgoing edges from u are considered.
	 */
	template<typename L> void forWeightedNeighborsOf(node u, L handle) const;

	/**
	 * Iterate over all outgoing edges of a node and call handler (lamdba closure).
	 * For undirected graphs all edges incident to u are also outgoing edges.
	 */
	template<typename L> void forEdgesOf(node u, L handle) const;

	/**
	 * Iterate over all outgoing edges of a node and call handler (lamdba closure).
	 * For undirected graphs all edges incident to u are also outgoing edges.
	 *
	 * Handle takes parameters (u, v, w) where w is the edge weight.
	 */
	template<typename L> void forWeightedEdgesOf(node u, L handle) const;

	/**
	 * Iterate over all neighbors of a node and call handler (lamdba closure).
	 * For directed graphs only incoming edges from u are considered.
	 */
	template<typename L> void forInNeighborsOf(node u, L handle) const;

	/**
	 * For directed graphs only incoming edges from u are considered.
	 * Iterate over all edge weights of a node and call handler (lamdba closure).
	 */
	template<typename L> void forWeightedInNeighborsOf(node u, L handle) const;

	/**
	 * Iterate over all incoming edges of a node and call handler (lamdba closure).
	 * For undirected graphs all edges incident to u are also incoming edges.
	 */
	template<typename L> void forInEdgesOf(node u, L handle) const;

	/**
	 * Iterate over all incoming edges of a node and call handler (lamdba closure).
	 * For undirected graphs all edges incident to u are also incoming edges.
	 *
	 * Handle takes parameters (u, v, w) where w is the edge weight.
	 */
	template<typename L> void forWeightedInEdgesOf(node u, L handle) const;


	/** REDUCTION ITERATORS **/
	
	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForNodes(L handle) const;

	/**
	 * Iterate in parallel over all edges and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForWeightedEdges(L handle) const;


	/** GRAPH SEARCHES **/

	template<typename L> void BFSfrom(node r, L handle) const;

	template<typename L> void BFSEdgesfrom(node r, L handle) const;

	template<typename L> void DFSfrom(node r, L handle) const;

	template<typename L> void DFSEdgesfrom(node r, L handle) const;
};

/** NODE ITERATORS **/

template<typename L>
void Graph::forNodes(L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
void Graph::parallelForNodes(L handle) const {
	#pragma omp parallel for
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename C, typename L>
void Graph::forNodesWhile(C condition, L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			if (!condition()) {
				break;
			}
			handle(v);
		}
	}
}

template<typename L>
void Graph::forNodesInRandomOrder(L handle) const {
	std::vector<node> randVec = nodes();
	random_shuffle(randVec.begin(), randVec.end());
	for (node v : randVec) {
		handle(v);
	}
}

template<typename L>
void Graph::balancedParallelForNodes(L handle) const {
	#pragma omp parallel for schedule(guided) // TODO: define min block size (and test it!)
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
void Graph::forNodePairs(L handle) const {
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

template<typename L>
void Graph::parallelForNodePairs(L handle) const {
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

template<typename L>
void Graph::forEdges(L handle) const {
	forWeightedEdges([&handle](node u, node v, edgeweight ew) { handle(u, v); });
}

template<typename L>
void Graph::parallelForEdges(L handle) const {
	parallelForWeightedEdges([&handle](node u, node v, edgeweight ew) { handle(u, v); });
}

template<typename L>
void Graph::forWeightedEdges(L handle) const {
	switch (weighted + 2 * directed) {
		case 0: // unweighted, undirected
			for (node u = 0; u < z; ++u) {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					node v = outEdges[u][i];
					// undirected, do not iterate over edges twice
					// {u, v} instead of (u, v); if v == none, u > v is not fulfilled
					if (u >= v) {
						edgeweight ew = defaultEdgeWeight;
						handle(u, v, ew);
					}
				}
			}
			break;
		
		case 1: // weighted,   undirected
			for (node u = 0; u < z; ++u) {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					node v = outEdges[u][i];
					// undirected, do not iterate over edges twice
					// {u, v} instead of (u, v); if v == none, u > v is not fulfilled
					if (u >= v) {
						edgeweight ew = outEdgeWeights[u][i];
						handle(u, v, ew);
					}
				}
			}
			break;

		case 2: // unweighted, directed
			for (node u = 0; u < z; ++u) {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					node v = outEdges[u][i];
					if (v != none) {
						edgeweight ew = defaultEdgeWeight;
						handle(u, v, ew);
					}
				}
			}
			break;

		case 3: // weighted,   directed
			for (node u = 0; u < z; ++u) {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					node v = outEdges[u][i];
					if (v != none) {
						edgeweight ew = outEdgeWeights[u][i];
						handle(u, v, ew);
					}
				}
			}
			break;
	}
}

template<typename L>
void Graph::parallelForWeightedEdges(L handle) const {
	switch (weighted + 2 * directed) {
		case 0: // unweighted, undirected
			#pragma omp parallel for
			for (node u = 0; u < z; ++u) {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					node v = outEdges[u][i];
					// undirected, do not iterate over edges twice
					// {u, v} instead of (u, v); if v == none, u > v is not fulfilled
					if (u >= v) {
						edgeweight ew = defaultEdgeWeight;
						handle(u, v, ew);
					}
				}
			}
			break;
		
		case 1: // weighted,   undirected
			#pragma omp parallel for
			for (node u = 0; u < z; ++u) {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					node v = outEdges[u][i];
					// undirected, do not iterate over edges twice
					// {u, v} instead of (u, v); if v == none, u > v is not fulfilled
					if (u >= v) {
						edgeweight ew = outEdgeWeights[u][i];
						handle(u, v, ew);
					}
				}
			}
			break;

		case 2: // unweighted, directed
			#pragma omp parallel for
			for (node u = 0; u < z; ++u) {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					node v = outEdges[u][i];
					if (v != none) {
						edgeweight ew = defaultEdgeWeight;
						handle(u, v, ew);
					}
				}
			}
			break;

		case 3: // weighted,   directed
			#pragma omp parallel for
			for (node u = 0; u < z; ++u) {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					node v = outEdges[u][i];
					if (v != none) {
						edgeweight ew = outEdgeWeights[u][i];
						handle(u, v, ew);
					}
				}
			}
			break;
	}
}

template<typename L>
void Graph::forEdgesWithAttribute_double(int attrId, L handle) const {
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < outEdges[u].size(); ++i) {
			node v = outEdges[u][i];
			if (directed) {
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

template<typename L>
void Graph::forNeighborsOf(node u, L handle) const {
	forWeightedEdgesOf(u, [&handle](node u, node v, edgeweight ew) { handle(v); });
}

template<typename L>
void Graph::forWeightedNeighborsOf(node u, L handle) const {
	forWeightedEdgesOf(u, [&handle](node u, node v, edgeweight ew) { handle(v, ew); });
}

template<typename L>
void Graph::forEdgesOf(node u, L handle) const {
	forWeightedEdgesOf(u, [&handle](node u, node v, edgeweight ew) { handle(u, v); });
}

template<typename L>
void Graph::forWeightedEdgesOf(node u, L handle) const {
	if (weighted) {
		for (index i = 0; i < outEdges[u].size(); i++) {
			node v = outEdges[u][i];
			if (v != none) {
				edgeweight ew = outEdgeWeights[u][i];
				handle(u, v, ew);
			}
		}	
	} else {
		for (index i = 0; i < outEdges[u].size(); i++) {
			node v = outEdges[u][i];
			if (v != none) {
				edgeweight ew = defaultEdgeWeight;
				handle(u, v, ew);
			}
		}
	}
}

template<typename L>
void Graph::forInNeighborsOf(node u, L handle) const {
	forWeightedInEdgesOf(u, [&handle](node u, node v, edgeweight ew) { handle(v); });
}

template<typename L>
void Graph::forWeightedInNeighborsOf(node u, L handle) const {
	forWeightedInEdgesOf(u, [&handle](node u, node v, edgeweight ew) { handle(v, ew); });
}

template<typename L>
void Graph::forInEdgesOf(node u, L handle) const {
	forWeightedInEdgesOf(u, [&handle](node u, node v, edgeweight ew) { handle(u, v); });
}

template<typename L>
void Graph::forWeightedInEdgesOf(node u, L handle) const {
	if (!directed) {
		forWeightedEdgesOf(u, handle);
		return;
	}
	if (weighted) {
		for (index i = 0; i < inEdges[u].size(); i++) {
			node v = inEdges[u][i];
			if (v != none) {
				edgeweight ew = inEdgeWeights[u][i];
				handle(u, v, ew);
			}
		}	
	} else {
		for (index i = 0; i < inEdges[u].size(); i++) {
			node v = inEdges[u][i];
			if (v != none) {
				edgeweight ew = defaultEdgeWeight;
				handle(u, v, ew);
			}
		}
	}
}

/** REDUCTION ITERATORS **/

template<typename L>
double Graph::parallelSumForNodes(L handle) const {
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			sum += handle(v);
		}
	}
	return sum;
}


/** GRAPH SEARCHES **/

template<typename L>
void Graph::BFSfrom(node r, L handle) const {
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

template<typename L>
void Graph::BFSEdgesfrom(node r, L handle) const {
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

template<typename L>
void Graph::DFSfrom(node r, L handle) const {
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

template<typename L>
void Graph::DFSEdgesfrom(node r, L handle) const {
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

} /* namespace NetworKit */

#endif /* BASICGRAPH_H_ */
