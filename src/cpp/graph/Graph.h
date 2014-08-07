/*
 * BasicGraph.h
 *
 *  Created on: 01.06.2014
 *      Author: Christian Staudt (christian.staudt@kit.edu), Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef GRAPH_H_
#define GRAPH_H_

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
#include "../auxiliary/Log.h"



namespace NetworKit {

	class ParallelPartitionCoarsening; // forward declaration for friend class
	class GraphBuilder; // forward declaration

/**
 * @ingroup graph
 * A graph (with optional weights) and parallel iterator methods.
 */
class Graph final {

	friend class ParallelPartitionCoarsening;
	friend class GraphBuilder;

private:
	// graph attributes
	count id; //!< unique graph id, starts at 0
	std::string name; //!< name of the graph, initially G#ID

	// scalars
	count n; //!< current number of nodes
	count m; //!< current number of edges
	node z; //!< current upper bound of node ids, z will be the id of the next node
	count t; //!< current time step

	bool weighted; //!< true if the graph is weighted, false otherwise
	bool directed; //!< true if the graph is directed, false otherwise

	// per node data
	std::vector<bool> exists; //!< exists[v] is true if node v has not been removed from the graph
	Coordinates<float> coordinates; //!< coordinates of nodes (if present)

	std::vector<count> inDeg; //!< only used for directed graphs, number of edges incoming per node
	std::vector<count> outDeg; //!< degree of every node, zero if node was removed. For directed graphs only outgoing edges count

	std::vector< std::vector<node> > inEdges; //!< only used for directed graphs, inEdges[v] contains all nodes u that have an edge (u, v)
	std::vector< std::vector<node> > outEdges; //!< (outgoing) edges, for each edge (u, v) v is saved in outEdges[u] and for undirected also u in outEdges[v]

	std::vector< std::vector<edgeweight> > inEdgeWeights; //!< only used for directed graphs, same schema as inEdges
	std::vector< std::vector<edgeweight> > outEdgeWeights; //!< same schema (and same order!) as outEdges

	// user-defined edge attributes
	// attribute maps storage
	std::vector<std::vector<std::vector<double> > > edgeMaps_double; //!< contains edge maps (u, v) -> double
	// default values
	std::vector<double> edgeAttrDefaults_double; //!< stores default value for edgeMaps_double[i] at index i

	/**
	 * Returns the next unique graph id.
	 */
	count getNextGraphId();

	/**
	 * Returns the index of node u in the array of incoming edges of node v. (for directed graphs inEdges is searched, while for indirected outEdges is searched, which gives the same result as indexInOutEdgeArray).
	 */
	index indexInInEdgeArray(node v, node u) const;

	/**
	 * Returns the index of node v in the array of outgoing edges of node u.
	 */
	index indexInOutEdgeArray(node u, node v) const;


public:

	/**
	 * Create a graph of @a n nodes. The graph has assignable edge weights if @a weighted is set to <code>true</code>.
	 * If @a weighted is set to <code>false</code> each edge has edge weight 1.0 and any other weight assignment will
	 * be ignored.
	 * @param n Number of nodes.
	 * @param weighted If set to <code>true</code>, the graph has edge weights.
	 * @param directed If set to @c true, the graph will be directed.
	 */
	Graph(count n = 0, bool weighted = false, bool directed = false);

	Graph(const Graph& G, bool weighted, bool directed);


	/**
	 * Create a graph as copy of @a other.
	 * @param other The graph to copy.
	 */
	Graph(const Graph& other) = default;

	/** Default move constructor */
	Graph(Graph&& other) = default;

	/** Default destructor */
	~Graph() = default;

	/** Default move assignment operator */
	Graph& operator=(Graph&& other) = default;

	/** Default copy assignment operator */
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
	 * Try to save some memory by shrinking internal data structures of the graph. Only run this
	 * once you finished editing the graph. Otherwise it will cause unnecessary reallocation of
	 * memory.
	 */
	void shrinkToFit();

	/**
	 * Set name of graph to @a name.
	 * @param name The name.
	 */
	void setName(std::string name) { this->name = name; }

	/*
	 * Returns the name of the graph.
	 * @return The name of the graph.
	 */
	std::string getName() const { return name; }


	/**
	 * Returns a string representation of the graph.
	 * @return A string representation.
	 */
	std::string toString() const;


	/* COPYING */

	/*
	* Copies all nodes to a new graph
	* @return graph with the same nodes.
	*/
	Graph copyNodes() const;


	/* NODE MODIFIERS */

	/**
	 * Add a new node to the graph and return it.
	 * @return The new node.
	 */
	node addNode();

	/**
	 * Add a new node to the graph with coordinates @a x and @y and return it.
	 */
	node addNode(float x, float y);

	/**
	 * Remove an isolated node @a v from the graph.
	 *
	 * @param u Node.
	 * @note Although it would be convenient to remove all incident edges at the same time,
	 * this causes complications for dynamic applications. Therefore, removeNode is an
	 * atomic event. All incident edges need to be removed first and an exception is thrown
	 * otherwise.
	 */
	void removeNode(node v);

	/**
	 * Check if node @a v exists in the graph.
	 *
	 * @param v Node.
	 * @return @c true if @a v exists, @c false otherwise.
	 */
	bool hasNode(node v) const { return (v < z) && this->exists[v];	}


	/** NODE PROPERTIES **/

	/**
	 * Returns the number of outgoing neighbors of @a v.
	 *
	 * @param v Node.
	 * @return The number of outgoing neighbors.
	 */
	count degree(node v) const { return outDeg[v]; }

	/**
	 * Get the number of incoming neighbors of @a v.
	 *
	 * @param v Node.
	 * @return The number of incoming neighbors.
	 * @note If the graph is not directed, the outgoing degree is returned.
	 */
	count degreeIn(node v) const { return directed ? inDeg[v] : outDeg[v]; }

	/**
	 * Get the number of outgoing neighbors of @a v.
	 *
	 * @param v Node.
	 * @return The number of outgoing neighbors.
	 */
	count degreeOut(node v) const { return outDeg[v]; }

	/**
	 * Check whether @a v is isolated, i.e. degree is 0.
	 * @param v Node.
	 * @return @c true if the node is isolated (= degree is 0)
	 */
	bool isIsolated(node v) const { return outDeg[v] == 0 && (!directed || inDeg[v] == 0); }


	/**
	 * Returns the weighted degree of @a v.
	 *
	 * @param v Node.
	 * @return Weighted degree of @a v.
	 * @note For directed graphs this is the sum of weights of all outgoing edges of @a v.
	 */
	edgeweight weightedDegree(node v) const;

	/**
	 * Returns the volume of the @a v, which is the weighted degree with self-loops counted twice.
	 *
	 * @param v Node.
	 * @return The volume of the @a v.
	 */
	edgeweight volume(node v) const;

	/**
	 * Returns a random node of the graph.
	 * @return A random node.
	 */
	node randomNode() const;

	/**
	 * Returns a random neighbor of @a u and @c none if degree is zero.
	 *
	 * @param u Node.
	 * @return A random neighbor of @a u.
	 */
	node randomNeighbor(node u) const;


	/* EDGE MODIFIERS */

	/**
	 * Insert an edge between the nodes @a u and @a v. If the graph is weighted you can optionally
	 * set a weight for this edge. The default weight is 1.0.
	 * @param u Endpoint of edge.
	 * @param v Endpoint of edge.
	 * @param weight Optional edge weight.
	 */
	void addEdge(node u, node v, edgeweight ew = defaultEdgeWeight);

	/**
	 * Removes the undirected edge {@a u,@a v}.
	 * @param u Endpoint of edge.
	 * @param v Endpoint of edge.
	 */
	void removeEdge(node u, node v);

	/**
	 * Checks if undirected edge {@a u,@a v} exists in the graph.
	 * @param u Endpoint of edge.
	 * @param v Endpoint of edge.
	 * @return <code>true</code> if the edge exists, <code>false</code> otherwise.
	 */
	bool hasEdge(node u, node v) const;

	/**
	 * Returns a random edge. By default a random node u is chosen and then some random neighbor v. So the probability of choosing (u, v) highly
	 * depends on the degree of u.
	 * Setting uniformDistribution to true, will give you a real uniform distributed edge, but will be very slow. So only use uniformDistribution
	 * for single calls outside of any loops.
	 */
	std::pair<node, node> randomEdge(bool uniformDistribution = false) const;

	/**
	 * Returns a vector with nr random edges. The edges are chosen uniform random.
	 */
	std::vector< std::pair<node, node> > randomEdges(count nr) const;

	/* GLOBAL PROPERTIES */

	/**
	 * Returns <code>true</code> if this graph supports edge weights other than 1.0.
	 * @return <code>true</code> if this graph supports edge weights other than 1.0.
	 */
	bool isWeighted() const { return weighted; }

	/**
	 * Return @c true if this graph supports directed edges.
	 * @return @c true if this graph supports directed edges.
	 */
	bool isDirected() const { return directed; }

	/**
	 * Return <code>true</code> if graph contains no nodes.
	 * @return <code>true</code> if graph contains no nodes.
	 */
	bool isEmpty() const { return n == 0; }

	/**
	 * Return the number of nodes in the graph.
	 * @return The number of nodes.
	 */
	count numberOfNodes() const { return n; }

	/**
	 * Return the number of edges in the graph.
	 * @return The number of edges.
	 */
	count numberOfEdges() const { return m; }

	/**
	 * Return the number of loops {v,v} in the graph.
	 * @return The number of loops.
	 * @note This involves calculation, so store result if needed multiple times.
	 */
	count numberOfSelfLoops() const;

 	/**
	 * Get an upper bound for the node ids in the graph.
	 * @return An upper bound for the node ids.
	 */
	index upperNodeIdBound() const { return z; }

	/**
	 * Check for invalid graph states, such as multiedges
	 */
	bool consistencyCheck() const;


	/* DYNAMICS */

	/**
	 * Trigger a time step - increments counter.
	 */
	void timeStep() { t++; }

	/**
	 * Get time step counter.
	 * @return Time step counter.
	 */
	count time() { return t; }


	/* COORDINATES */

	/**
	 * Sets the coordinate of @a v to @a value.
	 *
	 * @param v Node.
	 * @param value The coordinate of @a v.
	 */
	void setCoordinate(node v, Point<float> value) { coordinates.setCoordinate(v, value); }


	/**
	 * Get the coordinate of @a v.
	 * @param v Node.
	 * @return The coordinate of @a v.
	 */
	Point<float>& getCoordinate(node v) { return coordinates.getCoordinate(v); }

	/**
	 * Get minimum coordinate of all coordinates with respect to dimension @a dim.
	 * @param dim The dimension to search for minimum.
	 * @return The minimum coordinate in dimension @a dim.
	 */
	float minCoordinate(count dim) { return coordinates.minCoordinate(dim); }

	/**
	 * Get maximum coordinate of all coordinates with respect to dimension @a dim.
	 * @param dim The dimension to search for maximum.
	 * @return The maximum coordinate in dimension @a dim.
	 */
	float maxCoordinate(count dim) { return coordinates.maxCoordinate(dim); }

	/**
	 * Initializes the coordinates for the nodes in graph.
	 * @note This has to be called once and before you set coordinates. Call this method again if new nodes have
	 * been added.
	 */
	void initCoordinates() { coordinates.init(z); }


	/* EDGE ATTRIBUTES */

	/**
	 * Return edge weight of edge {@a u,@a v}. Returns 0 if edge does not exist.
	 *
	 * @param u Endpoint of edge.
	 * @param v Endpoint of edge.
	 * @return Edge weight of edge {@a u,@a v} or 0 if edge does not exist.
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
	 * Add new edge map for an attribute of type <code>double</code> with @a defaultValue.
	 *
	 * @param defaultValue The default value if no other value for this attribute is specified.
	 * @return The attribute id of this map.
	 */
	int addEdgeAttribute_double(double defaultValue);


	/**
	 * Returns the attribute of type <code>double</code> with @a attrId for edge {@a u,@a v}.
	 *
	 * @param[in]	u	Endpoint of edge.
	 * @param[in]	v	Endpoint of edge.
	 * @param[in]	attrId	Attribute id.
	 * @return Attribute with @a attrId for edge {@a u,@a v}.
	 */
	double attribute_double(node u, node v, int attrId) const;


	/**
	 * Set edge attribute @a attr of type <code>double</code> with @a attrId of edge {@a u,@a v}. If the edge
	 * does not exist, it will be inserted.
	 *
	 * @param[in]	u	Endpoint of edge.
	 * @param[in]	v	Endpoint of edge.
	 * @param[in]	attrId Attribute id.
	 * @param[in]	attr	Edge attribute.
	 */
	void setAttribute_double(node u, node v, int attrId, double attr);


	/* SUMS */

	/**
	 * Returns the sum of all edge weights.
	 * @return The sum of all edge weights.
	 */
	edgeweight totalEdgeWeight() const;


	/* Collections */

	/**
	 * Get list of all nodes.
	 * @return List of all nodes.
	 */
	std::vector<node> nodes() const;

	/**
	 * Get list of edges as node pairs.
	 * @return List of edges as node pairs.
	 */
	std::vector<std::pair<node, node> > edges() const;

	/**
	 * Get list of neighbors of @a u.
	 *
	 * @param u Node.
	 * @return List of neighbors of @a u.
	 */
	std::vector<node> neighbors(node u) const;


	/* NODE ITERATORS */

	/**
	 * Iterate over all nodes of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void forNodes(L handle) const;

	/**
	 * Iterate randomly over all nodes of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void parallelForNodes(L handle) const;

	/** Iterate over all nodes of the graph and call @a handle (lambda closure) as long as @a condition remains true.
	 * This allows for breaking from a node loop.
	 *
	 * @param condition Returning <code>false</code> breaks the loop.
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename C, typename L> void forNodesWhile(C condition, L handle) const;

	/**
	 * Iterate randomly over all nodes of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void forNodesInRandomOrder(L handle) const;

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 * Using schedule(guided) to remedy load-imbalances due to e.g. unequal degree distribution.
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void balancedParallelForNodes(L handle) const;


	/**
	 * Iterate over all undirected pairs of nodes and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void forNodePairs(L handle) const;


	/**
	 * Iterate over all undirected pairs of nodes in parallel and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void parallelForNodePairs(L handle) const;


	/* EDGE ITERATORS */

	/**
	 * Iterate over all edges of the const graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void forEdges(L handle) const;

	/**
	 * Iterate in parallel over all edges of the const graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void parallelForEdges(L handle) const;

	/**
	 * Iterate over all edges of the const graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node, edgeweight)</code>.
	 */
	template<typename L> void forWeightedEdges(L handle) const;

	/**
	 * Iterate over all edges of the const graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node, edgeweight)</code>.
	 */
	template<typename L> void parallelForWeightedEdges(L handle) const;


	/**
	 * Iterate over all edges of the const graph and call @a handle (lambda closure).
	 *
	 *	@param[in] 	attrId		Attribute id.
	 *  @param[in] 	handle 		Takes parameters <code>(node, node, double)</code>.
	 */
	template<typename L> void forEdgesWithAttribute_double(int attrId, L handle) const;

	/* NEIGHBORHOOD ITERATORS */

	/**
	 * Iterate over all neighbors of a node and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameter <code>(node)</code> which is a neighbor of @a u.
	 * @note For directed graphs only outgoing edges from @a u are considered.
	 * A node is its own neighbor if there is a self-loop.
	 *
	 */
	template<typename L> void forNeighborsOf(node u, L handle) const;

	/**
	 * Iterate over all edge weights of a node and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameters <code>(node, const edgeweight)</code> where node is a neighbor of @a u.
	 * @note For directed graphs only outgoing edges from u are considered.
	 * A node is its own neighbor if there is a self-loop.
	 */
	template<typename L> void forWeightedNeighborsOf(node u, L handle) const;


	/**
	 * Iterate over all incident edges of a node and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameters <code>(node, node)</code> where the first node is @a u and the second is a neighbor of @a u.
	 * @note For undirected graphs all edges incident to @a u are also outgoing edges.
	 */
	template<typename L> void forEdgesOf(node u, L handle) const;

	/**
	 * Iterate over all outgoing edges of a node and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameters <code>(node, node, edgeweight)</code> where the first node is @a u and the second is
	 * a neighbor of @a u.
	 * @note For undirected graphs all edges incident to @a u are also outgoing edges.
	 */
	template<typename L> void forWeightedEdgesOf(node u, L handle) const;

	/**
	 * Iterate over all neighbors of a node and call handler (lamdba closure).
	 * For directed graphs only incoming edges from u are considered.
	 */
	template<typename L> void forInNeighborsOf(node u, L handle) const;

	/**
	 * Iterate over all edge weights of a node and call handler (lamdba closure).
	 *
	 * @note For directed graphs only incoming edges from u are considered.
	 */
	template<typename L> void forWeightedInNeighborsOf(node u, L handle) const;

	/**
	 * Iterate over all incoming edges of a node and call handler (lamdba closure).
	 * @note For undirected graphs all edges incident to u are also incoming edges.
	 */
	template<typename L> void forInEdgesOf(node u, L handle) const;

	/**
	 * Iterate over all incoming edges of a node and call handler (lamdba closure).
	 * @note For undirected graphs all edges incident to u are also incoming edges.
	 *
	 * Handle takes parameters (u, v, w) where w is the edge weight.
	 */
	template<typename L> void forWeightedInEdgesOf(node u, L handle) const;


	/* REDUCTION ITERATORS */

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForNodes(L handle) const;

	/**
	 * Iterate in parallel over all edges and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForEdges(L handle) const;

	/**
	 * Iterate in parallel over all edges and sum (reduce +) the values returned by the handler
	 *
	 * @param handle Takes parameter <code>(node)</code> and returns <code>double</code>.
	 */
	template<typename L> double parallelSumForWeightedEdges(L handle) const;


	/* GRAPH SEARCHES */

	/**
	 * Iterate over nodes in breadth-first search order starting from r until connected component
	 * of r has been visited.
	 *
	 * @param r Node.
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void BFSfrom(node r, L handle) const;
	template<typename L> void BFSfrom(std::vector<node> &startNodes, L handle) const;


	template<typename L> void BFSEdgesFrom(node r, L handle) const;

	/**
	 * Iterate over nodes in depth-first search order starting from r until connected component
	 * of r has been visited.
	 *
	 * @param r Node.
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void DFSfrom(node r, L handle) const;


	template<typename L> void DFSEdgesFrom(node r, L handle) const;


	/* SPECIAL */


	/**
	* Treat the adjacency datastructure as an undirected graph.
	*/
	void treatAsUndirected();
};

/* NODE ITERATORS */

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
	#pragma omp parallel for schedule(guided)
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


/* EDGE ITERATORS */

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


/* NEIGHBORHOOD ITERATORS */

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
	forWeightedInEdgesOf(u, [&handle](node u, node v, edgeweight ew) { handle(u); });
}

template<typename L>
void Graph::forWeightedInNeighborsOf(node u, L handle) const {
	forWeightedInEdgesOf(u, [&handle](node u, node v, edgeweight ew) { handle(u, ew); });
}

template<typename L>
void Graph::forInEdgesOf(node u, L handle) const {
	forWeightedInEdgesOf(u, [&handle](node u, node v, edgeweight ew) { handle(u, v); });
}

template<typename L>
void Graph::forWeightedInEdgesOf(node u, L handle) const {
	if (directed) {
		if (weighted) {
			for (index i = 0; i < inEdges[u].size(); i++) {
				node v = inEdges[u][i];
				if (v != none) {
					edgeweight ew = inEdgeWeights[u][i];
					handle(v, u, ew);
				}
			}
		} else {
			for (index i = 0; i < inEdges[u].size(); i++) {
				node v = inEdges[u][i];
				if (v != none) {
					edgeweight ew = defaultEdgeWeight;
					handle(v, u, ew);
				}
			}
		}
	} else {
		// use outEdges
		// the following dose the same as calling forWeightedEdgesOf(u, [&handle](node u, node v, edgeweight ew) { handle(v, u, ew); });
		// but I think it is better to write it out here
		if (weighted) {
			for (index i = 0; i < outEdges[u].size(); i++) {
				node v = outEdges[u][i];
				if (v != none) {
					edgeweight ew = outEdgeWeights[u][i];
					handle(v, u, ew);
				}
			}
		} else {
			for (index i = 0; i < outEdges[u].size(); i++) {
				node v = outEdges[u][i];
				if (v != none) {
					edgeweight ew = defaultEdgeWeight;
					handle(v, u, ew);
				}
			}
		}
	}
}

/* REDUCTION ITERATORS */

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

template<typename L>
double Graph::parallelSumForEdges(L handle) const {
	return parallelSumForWeightedEdges([&handle](node u, node v, edgeweight ew) { return handle(u, v); });
}

template<typename L>
double Graph::parallelSumForWeightedEdges(L handle) const {
	double sum = 0.0;
	switch (weighted + 2 * directed) {
		case 0: // unweighted, undirected
			#pragma omp parallel for reduction(+:sum)
			for (node u = 0; u < z; ++u) {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					node v = outEdges[u][i];
					// undirected, do not iterate over edges twice
					// {u, v} instead of (u, v); if v == none, u > v is not fulfilled
					if (u >= v) {
						edgeweight ew = defaultEdgeWeight;
						sum += handle(u, v, ew);
					}
				}
			}
			break;

		case 1: // weighted,   undirected
			#pragma omp parallel for reduction(+:sum)
			for (node u = 0; u < z; ++u) {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					node v = outEdges[u][i];
					// undirected, do not iterate over edges twice
					// {u, v} instead of (u, v); if v == none, u > v is not fulfilled
					if (u >= v) {
						edgeweight ew = outEdgeWeights[u][i];
						sum +=handle(u, v, ew);
					}
				}
			}
			break;

		case 2: // unweighted, directed
			#pragma omp parallel for reduction(+:sum)
			for (node u = 0; u < z; ++u) {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					node v = outEdges[u][i];
					if (v != none) {
						edgeweight ew = defaultEdgeWeight;
						sum +=handle(u, v, ew);
					}
				}
			}
			break;

		case 3: // weighted,   directed
			#pragma omp parallel for reduction(+:sum)
			for (node u = 0; u < z; ++u) {
				for (index i = 0; i < outEdges[u].size(); ++i) {
					node v = outEdges[u][i];
					if (v != none) {
						edgeweight ew = outEdgeWeights[u][i];
						sum +=handle(u, v, ew);
					}
				}
			}
			break;
	}
	return sum;
}


/* GRAPH SEARCHES */

template<typename L>
void Graph::BFSfrom(node r, L handle) const {
	std::vector<node> startNodes(1, r);
	BFSfrom(startNodes, handle);
}

template<typename L>
void Graph::BFSfrom(std::vector<node> &startNodes, L handle) const {
	std::vector<bool> marked(z);
	std::queue<node> q, qNext;
	count dist = 0;
	// enqueue start nodes
	for (node u : startNodes) {
		q.push(u);
		marked[u] = true;
	}
	do {
		node u = q.front();
		q.pop();
		// apply function
		handle(u, dist);
		forNeighborsOf(u, [&](node v) {
			if (!marked[v]) {
				qNext.push(v);
				marked[v] = true;
			}
		});
		if (q.empty() && !qNext.empty()) {
			q.swap(qNext);
			++dist;
		}
	} while (!q.empty());
}

template<typename L>
void Graph::BFSEdgesFrom(node r, L handle) const {
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
void Graph::DFSEdgesFrom(node r, L handle) const {
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

#endif /* GRAPH_H_ */
