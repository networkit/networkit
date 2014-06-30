/*
 * Graph.h
 *
 *  Created on: 04.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu), Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <functional>
#include <cassert>
#include <vector>
#include <cinttypes>
#include <string>
#include <queue>
#include <stack>
#include <stdexcept>
#include <map>
#include <set>
#include <sstream>
#include <limits>
#include <cstdint>
#include <algorithm>
// #include <tbb/concurrent_vector.h>

#include "Coordinates.h"
#include "../auxiliary/Log.h"
#include "../Globals.h"
#include "../viz/Point.h"


namespace NetworKit {

/** Typedefs **/

typedef uint64_t index; // more expressive name for an index into an array
typedef uint64_t count; // more expressive name for an integer quantity
typedef index node; // node indices are 0-based
typedef double edgeweight; // edge weight type

constexpr index none = std::numeric_limits<index>::max();

//#define StdVector std::vector // TODO: test tbb::concurrent_vector
template <typename T> using StdVector = std::vector<T>;

/**
 * @ingroup graph
 * An undirected graph (with optional weights) and parallel iterator methods.
 */
class Graph final {

protected:

	// scalars
	count n; //!< current number of nodes
	count m; //!< current number of edges
	node z; //!< current upper bound of node ids
	count t; //!< current time step
	bool weighted; //!< true if this graph supports edge weights other than 1.0

	// per node data
	StdVector<count> deg; //!< degree of each node (size of neighborhood)
	StdVector<bool> exists; //!< exists[v] is true if node v has not been removed from the graph
	Coordinates<float> coordinates; //!< coordinates of nodes (if present)

	// per edge data
	StdVector<StdVector<node> > adja; //!< neighbors/adjacencies
	StdVector<StdVector<edgeweight> > eweights; //!< edge weights

	// graph attributes
	std::string name;

	// user-defined edge attributes

	//	attribute maps storage

	std::vector<std::vector<std::vector<double> > > edgeMaps_double; // contains edge maps (u, v) -> double

	// defaults

	std::vector<double> edgeAttrDefaults_double; // stores default value for edgeMaps_double[i] at index i


	/**
	 * Return the index of v in the adjacency array of u.
	 */
	index find(node u, node v) const;

public:

	// defaults
	static constexpr double defaultEdgeWeight = 1.00;
	static constexpr edgeweight nullWeight = 0.0;

	/** ATTRIBUTE ABSTRACT BASE CLASSES **/

	class NodeAttribute {
		// abstract
	};

	class EdgeAttribute {
		// abstract
	};


	/** GRAPH INTERFACE **/

	/**
	 * Create a graph of @a n nodes. The graph has assignable edge weights if @a weighted is set to <code>true</code>.
	 * If @a weighted is set to <code>false</code> each edge has edge weight 1.0 and any other weight assignment will
	 * be ignored.
	 * @param n Number of nodes.
	 * @param weighted If set to <code>true</code>, the graph has edge weights.
	 */
	Graph(count n=0, bool weighted=false);

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



	/**
	 * Set name of graph to @a name.
	 * @param name The name.
	 */
	void setName(std::string name);

	/*
	 * Returns the name of the graph.
	 * @return The name of the graph.
	 */
	std::string getName();

	/**
	 * Returns <code>true</code> if this graph supports edge weights other than 1.0.
	 * @return <code>true</code> if this graph supports edge weights other than 1.0.
	 */
	bool isWeighted() const;

	/**
	 * Returns a string representation of the graph.
	 * @return A string representation.
	 */
	std::string toString();

	/**
	 * Insert an undirected edge between the nodes @a u and @a v. If the graph is weighted you can optionally
	 * set a weight for this edge. The default weight is 1.0.
	 * @param u Endpoint of edge.
	 * @param v Endpoint of edge.
	 * @param weight Optional edge weight.
	 */
	void addEdge(node u, node v, edgeweight weight = defaultEdgeWeight);

	/**
	 * Checks if undirected edge {@a u,@a v} exists in the graph.
	 * @param u Endpoint of edge.
	 * @param v Endpoint of edge.
	 * @return <code>true</code> if the edge exists, <code>false</code> otherwise.
	 */
	bool hasEdge(node u, node v) const;

	/**
	 * Removes the undirected edge {@a u,@a v}.
	 * @param u Endpoint of edge.
	 * @param v Endpoint of edge.
	 */
	void removeEdge(node u, node v);

	/**
	 * Merges edge {@a u,@a v} to become a supernode. Edges to u and v are
	 * rewired, multiple edges merged and their weights added.
	 * The vertex weights of @a u and @a v are added.
	 * A self-loop is only created if @a discardSelfLoop is set to <code>false</code>.
	 *
	 * @param u Endpoint of edge.
	 * @param v Endpoint of edge.
	 * @param discardSelfLoop If set to <code>true</code> (default) no self-loop is created.
	 * @return New node that has been created if @a u != @a v. Otherwise none.
	 */
	node mergeEdge(node u, node v, bool discardSelfLoop = true);

	/**
	 * Returns the number of neighbors of @a v.
	 *
	 * @param v Node.
	 * @return The number of neighbors.
	 */
	count degree(node v) const;

	/**
	 * Returns the smallest neighborhood size (does not have to be unique).
	 * @return The smallest neighborhood size.
	 */
	count minDegree() const;

	/**
	 * Returns the index of a node with smallest neighborhood size (does not have to be unique).
	 * @return The index of a node with smallest neighborhood size.
	 */
	index argminDegree() const;

	/**
	 * Returns the largest neighborhood size (does not have to be unique).
	 * @return The largest neighborhood size.
	 */
	count maxDegree() const;

	/**
	 * Returns the index of a node with largest neighborhood size (does not have to be unique).
	 * @return The index of a node with largest neighborhood size.
	 */
	index argmaxDegree() const;

	/**
	 * Returns the weighted degree of @a v.
	 *
	 * @param v Node.
	 * @return Weighted degree of @a v.
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
	 * Returns a random neighbor of @a v and <code>none</code> if degree is zero.
	 *
	 * @param v Node.
	 * @return A random neighbor of @a v.
	 */
	node randomNeighbor(node v) const;


    /**
     * Returns a random edge of the graph.
     * @return Random random edge.
     * @note Fast, but not uniformly random.
     */
    std::pair<node, node> randomEdge() const;	// TODO: implement uniformly random edge choice






	/** EDGE ATTRIBUTE GETTERS **/

	/**
	 * Return edge weight of edge {@a u,@a v}. Returns 0 if edge does not exist.
	 *
	 * @param u Endpoint of edge.
	 * @param v Endpoint of edge.
	 * @return Edge weight of edge {@a u,@a v} or 0 if edge does not exist.
	 */
	edgeweight weight(node u, node v) const;

	/**
	 * Returns the attribute of type <code>double</code> with @a attrId for edge {@a u,@a v}.
	 *
	 * @param[in]	u	Endpoint of edge.
	 * @param[in]	v	Endpoint of edge.
	 * @param[in]	attrId	Attribute id.
	 * @return Attribute with @a attrId for edge {@a u,@a v}.
	 */
	double attribute_double(node u, node v, int attrId) const;

	/**  EDGE ATTRIBUTE SETTERS */

	/**
	 * Set the weight of edge {@a u,@a v} to @a w. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	Endpoint of edge.
	 * @param[in]	v	Endpoint of edge.
	 * @param[in]	weight	Edge weight.
	 */
	void setWeight(node u, node v, edgeweight w);


	/**
	 * Increase the weight of edge {@a u,@a v} by @a w. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	Endpoint of edge.
	 * @param[in]	v	Endpoint of edge.
	 * @param[in]	weight	Edge weight.
	 */
	void increaseWeight(node u, node v, edgeweight w);

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

	/** SUMS **/

	/**
	 * Returns the sum of all edge weights.
	 * @return The sum of all edge weights.
	 */
	edgeweight totalEdgeWeight() const;


	/** NODE MODIFIERS **/

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
	 * Remove the isolated node @a u from the graph.
	 *
	 * @param u Node.
	 * @note Although it would be convenient to remove all incident edges at the same time,
	 * this causes complications for dynamic applications. Therefore, removeNode is an
	 * atomic event. All incident edges need to be removed first and an exception is thrown
	 * otherwise.
	 */
	void removeNode(node u);

	/**
	 * Check if node @a u exists in the graph.
	 *
	 * @param u Node.
	 * @return <code>true</code> if @a u exists, <code>false</code> otherwise.
	 */
	bool hasNode(node u) const;


	/** GLOBAL PROPERTIES **/

	/**
	 * Return <code>true</code> if graph contains no nodes.
	 * @return <code>true</code> if graph contains no nodes.
	 */
	bool isEmpty();

	/**
	 * Return the number of nodes in the graph.
	 * @return The number of nodes.
	 */
	count numberOfNodes() const;

	/**
	 * Return the number of edges in the graph.
	 * @return The number of edges.
	*/
	count numberOfEdges() const;

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
	index upperNodeIdBound() const;

	/** DYNAMICS **/

	/**
	 * Trigger a time step - increments counter.
	 */
	void timeStep();

	/**
	 * Get time step counter.
	 * @return Time step counter.
	 */
	count time();


	/** COORDINATES **/

	/**
	 * Sets the coordinate of @a v to @a value.
	 *
	 * @param v Node.
	 * @param value The coordinate of @a v.
	 */
	void setCoordinate(node v, Point<float> value) {
		coordinates.setCoordinate(v, value);
	}

	/**
	 * Get the coordinate of @a v.
	 * @param v Node.
	 * @return The coordinate of @a v.
	 */
	Point<float>& getCoordinate(node v) {
		return coordinates.getCoordinate(v);
	}

	/**
	 * Get minimum coordinate of all coordinates with respect to dimension @a dim.
	 * @param dim The dimension to search for minimum.
	 * @return The minimum coordinate in dimension @a dim.
	 */
	float minCoordinate(count dim) {
		return coordinates.minCoordinate(dim);
	}

	/**
	 * Get maximum coordinate of all coordinates with respect to dimension @a dim.
	 * @param dim The dimension to search for maximum.
	 * @return The maximum coordinate in dimension @a dim.
	 */
	float maxCoordinate(count dim) {
		return coordinates.maxCoordinate(dim);
	}

	/**
	 * Initializes the coordinates for the nodes in graph.
	 * @note This has to be called once and before you set coordinates. Call this method again if new nodes have
	 * been added.
	 */
	void initCoordinates() {
		coordinates.init(z);
	}

	/* ATTRIBUTES */

	/**
	 * Add new edge map for an attribute of type <code>double</code> with @a defaultValue.
	 *
	 * @param defaultValue The default value if no other value for this attribute is specified.
	 * @return The attribute id of this map.
	 */
	int addEdgeAttribute_double(double defaultValue);

	/* NODE ITERATORS */

	/**
	 * Iterate over all nodes of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void forNodes(L handle);

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
	template<typename L> void forNodesInRandomOrder(L handle);

	/**
	 * Iterate randomly over all nodes of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void forNodesInRandomOrder(L handle) const;

	/**
	 * Iterate over all nodes of the graph and call @a handle (lambda closure) as long as @a condition remains true.
	 * This allows for breaking from a node loop.
	 *
	 * @param condition Returning <code>false</code> breaks the loop.
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename C, typename L> void forNodesWhile(C condition, L handle);

	/**
	 * Iterate over all nodes of the graph and call @a handle (lambda closure) as long as @a condition remains true.
	 * This allows for breaking from a node loop.
	 *
	 * @param condition Returning <code>false</code> breaks the loop.
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename C, typename L> void forNodes(C condition, L handle) const;

	/**
	 * Iterate in parallel over all nodes of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void parallelForNodes(L handle);

	/**
	 * Iterate in parallel over all nodes of the graph and call @a handle (lambda closure).
	 * Using schedule(guided) to remedy load-imbalances due to e.g. unequal degree distribution.
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void balancedParallelForNodes(L handle);

	/**
	 * Iterate in parallel over all nodes of the graph and call @a handle (lambda closure).
	 * Using schedule(guided) to remedy load-imbalances due to e.g. unequal degree distribution.
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void balancedParallelForNodes(L handle) const;

	/**
	 * Iterate in parallel over all nodes of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void parallelForNodes(L handle) const;

	/**
	 * Iterate over all undirected pairs of nodes and call @a handle (lambda closure).
	 *
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void forNodePairs(L handle);

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
	template<typename L> void parallelForNodePairs(L handle);

	/**
	 * Iterate over all undirected pairs of nodes in parallel and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void parallelForNodePairs(L handle) const;

	/**
	 * Iterate over nodes in breadth-first search order starting from r until connected component
	 * of r has been visited.
	 *
	 * @param r Node.
	 * @param marked %Vector of node size initialized to 0.
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void breadthFirstNodesFrom(node r,
			std::vector<int>& marked, L handle);

	/**
	 * Iterate over nodes in breadth-first search order starting from r until connected component
	 * of r has been visited.
	 *
	 * @param r Node.
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void BFSfrom(node r, L handle);

	/**
	 * Iterate over nodes in depth-first search order starting from r until connected component
	 * of r has been visited.
	 *
	 * @param r Node.
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void DFSfrom(node r, L handle);

	/**
	 * Iterate over edges in breadth-first search order starting from node r until connected component
	 * of r has been visited.
	 *
	 * @note Not working yet.
	 */
	template<typename L> void breadthFirstEdgesFrom(node r, L handle);

	/**
	 * Iterate over all nodes of the graph and call @a handle (lambda closure).
	 *
	 * @param[in]	attrKey		attribute key
	 * @param[in]	handle		takes parameters (v, a) where a is a node attribute
	 *
	 * @note Not working yet.
	 */
	template<typename L> void forNodesWithAttribute(std::string attrKey,
			L handle);

	/* EDGE ITERATORS */

	/**
	 * Iterate over all edges of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void forEdges(L handle);

	/**
	 * Iterate over all edges of the const graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void forEdges(L handle) const;

	/**
	 * Iterate in parallel over all edges of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void parallelForEdges(L handle);

	/**
	 * Iterate in parallel over all edges of the const graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void parallelForEdges(L handle) const;

	/**
	 * Iterate over all edges of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node, edgeweight)</code>.
	 */
	template<typename L> void forWeightedEdges(L handle);

	/**
	 * Iterate over all edges of the const graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node, edgeweight)</code>.
	 */
	template<typename L> void forWeightedEdges(L handle) const;

	/**
	 * Iterate over all edges of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node, edgeweight)</code>.
	 */
	template<typename L> void parallelForWeightedEdges(L handle);

	/**
	 * Iterate over all edges of the const graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node, edgeweight)</code>.
	 */
	template<typename L> void parallelForWeightedEdges(L handle) const;

	/**
	 * Iterate over all edges of the graph and call @a handle (lambda closure).
	 *
	 *	@param[in] 	attrId		Attribute id.
	 *  @param[in] 	handle 		Takes parameters <code>(node, node, double)</code>.
	 */
	template<typename L> void forEdgesWithAttribute_double(int attrId,
			L handle);

	/**
	 * Iterate over all edges of the const graph and call @a handle (lambda closure).
	 *
	 *	@param[in] 	attrId		Attribute id.
	 *  @param[in] 	handle 		Takes parameters <code>(node, node, double)</code>.
	 */
	template<typename L> void forEdgesWithAttribute_double(int attrId,
			L handle) const;

	/* NEIGHBORHOOD ITERATORS */

	/**
	 * Iterate over all neighbors of a node and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameter <code>(node)</code> which is a neighbor of @a u.
	 * @note A node is its own neighbor if there is a self-loop.
	 */
	template<typename L> void forNeighborsOf(node u, L handle);

	/**
	 * Iterate over all neighbors of a node and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameter <code>(node)</code> which is a neighbor of @a u.
	 * @note A node is its own neighbor if there is a self-loop.
	 */
	template<typename L> void forNeighborsOf(node u, L handle) const;

	/**
	 * Iterate over all edge weights of a node and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameters <code>(node, const edgeweight)</code> where node is a neighbor of @a u.
	 * @note A node is its own neighbor if there is a self-loop.
	 */
	template<typename L> void forWeightedNeighborsOf(node u, L handle);

	/**
	 * Iterate over all edge weights of a node and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameters <code>(node, const edgeweight)</code> where node is a neighbor of @a u.
	 * @note A node is its own neighbor if there is a self-loop.
	 */
	template<typename L> void forWeightedNeighborsOf(node u, L handle) const;

	/**
	 * Iterate over all incident edges of a node and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameters <code>(node, node)</code> where the first node is @a u and the second is a neighbor of @a u.
	 */
	template<typename L> void forEdgesOf(node u, L handle);

	/**
	 * Iterate over all incident edges of a node and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameters <code>(node, node)</code> where the first node is @a u and the second is a neighbor of @a u.
	 */
	template<typename L> void forEdgesOf(node u, L handle) const;

	/**
	 * Iterate over all incident edges of a node in neighborhood-size-increasing order and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameters <code>(node, node)</code> where the first node is @a u and the second is a neighbor of @a u.
	 */
	template<typename L> void forEdgesOfInDegreeIncreasingOrder(node u, L handle);

	/**
	 * Iterate over all incident edges of a node in neighborhood-size-increasing order and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameters <code>(node, node)</code> where the first node is @a u and the second is a neighbor of @a u.
	 */
	template<typename L> void forEdgesOfInDegreeIncreasingOrder(node u, L handle) const;

	/**
	 * Iterate over all incident edges of a node and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameters <code>(node, node, edgeweight)</code> where the first node is @a u and the second is
	 * a neighbor of @a u.
	 */
	template<typename L> void forWeightedEdgesOf(node u, L handle);

	/**
	 * Iterate over all incident edges of a node and call @a handle (lamdba closure).
	 *
	 * @param u Node.
	 * @param handle Takes parameters <code>(node, node, edgeweight)</code> where the first node is @a u and the second is
	 * a neighbor of @a u.
	 */
	template<typename L> void forWeightedEdgesOf(node u, L handle) const;

	/** REDUCTION ITERATORS **/

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by @a handle
	 *
	 * @param handle Takes parameter <code>(node)</code> and returns <code>double</code>.
	 */
	template<typename L> double parallelSumForNodes(L handle);

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForNodes(L handle) const;

	/**
	 * Iterate in parallel over all edges and sum (reduce +) the values returned by the handler
	 *
	 * @param handle Takes parameter <code>(node)</code> and returns <code>double</code>.
	 */
	template<typename L> double parallelSumForWeightedEdges(L handle) const;


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
	std::vector<std::pair<node, node> > edges();


	/**
	 * Get list of neighbors of @a u.
	 *
	 * @param u Node.
	 * @return List of neighbors of @a u.
	 */
	std::vector<node> neighbors(node u);


};

} /* namespace NetworKit */

template<typename L>
inline void NetworKit::Graph::forNeighborsOf(node u, L handle) {
	for (node v : this->adja[u]) {
		if (v != none) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forNeighborsOf(node u, L handle) const {
	for (node v : this->adja[u]) {
		if (v != none) {
			handle(v);
		}
	}
}


template<typename L>
inline void NetworKit::Graph::forWeightedNeighborsOf(node u, L handle) {
	if (weighted) {
		for (index i = 0; i < (index) adja[u].size(); ++i) {
			node v = adja[u][i];
			if (v != none) {
				edgeweight ew = eweights[u][i];
				handle(v, ew);
				assert(ew == weight(u, v));
			}
		}
	} else {
		for (index i = 0; i < (index) adja[u].size(); ++i) {
			node v = adja[u][i];
			if (v != none) {
				handle(v, defaultEdgeWeight);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forWeightedNeighborsOf(node u, L handle) const {
	if (weighted) {
		for (index i = 0; i < (index) adja[u].size(); ++i) {
			node v = adja[u][i];
			if (v != none) {
				edgeweight ew = eweights[u][i];
				handle(v, ew);
				assert(ew == weight(u, v));
			}
		}
	} else {
		for (index i = 0; i < (index) adja[u].size(); ++i) {
			node v = adja[u][i];
			if (v != none) {
				handle(v, defaultEdgeWeight);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forNodes(L handle) {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forNodes(L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::parallelForNodes(L handle) {
#pragma omp parallel for
	for (node v = 0; v < z; ++v) {
		// call here
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::parallelForNodes(L handle) const {
#pragma omp parallel for
	for (node v = 0; v < z; ++v) {
		// call here
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::balancedParallelForNodes(L handle) {
#pragma omp parallel for schedule(guided) // TODO: define min block size (and test it!)
	for (node v = 0; v < z; ++v) {
		// call here
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::balancedParallelForNodes(L handle) const {
#pragma omp parallel for schedule(guided) // TODO: define min block size (and test it!)
	for (node v = 0; v < z; ++v) {
		// call here
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline double NetworKit::Graph::parallelSumForNodes(L handle) {
	double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (node v = 0; v < z; ++v) {
		// call here
		if (exists[v]) {
			sum += handle(v);
		}
	}
	return sum;
}

template<typename L>
inline double NetworKit::Graph::parallelSumForNodes(L handle) const {
	double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (node v = 0; v < z; ++v) {
		// call here
		if (exists[v]) {
			sum += handle(v);
		}
	}
	return sum;
}

template<typename L>
double NetworKit::Graph::parallelSumForWeightedEdges(L handle) const {
	double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < this->adja[u].size(); ++i) {
			node v = this->adja[u][i];
			edgeweight ew = this->eweights[u][i];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				sum += handle(u, v, ew);
			}
		}
	}
	return sum;
}

template<typename L>
inline void NetworKit::Graph::forEdges(L handle) {
	for (node u = 0; u < z; ++u) {
		for (node v : this->adja[u]) {
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forEdges(L handle) const {
	for (node u = 0; u < z; ++u) {
		for (node v : this->adja[u]) {
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::parallelForEdges(L handle) {
#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (node v : this->adja[u]) {
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::parallelForEdges(L handle) const {
#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (node v : this->adja[u]) {
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				handle(u, v);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forNodePairs(L handle) {
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				// call node pair function
				if (exists[v]) {
					handle(u, v);
				}
			}
		}

	}
}

template<typename L>
inline void NetworKit::Graph::forNodePairs(L handle) const {
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				// call node pair function
				if (exists[v]) {
					handle(u, v);
				}
			}
		}

	}
}

template<typename L>
inline void NetworKit::Graph::breadthFirstNodesFrom(node r,
		std::vector<int>& marked, L handle) {
	std::queue<node> q;
	q.push(r); // enqueue root
	marked[r] = 1;
	do {
		node u = q.front();
		q.pop();
		// apply function
		handle(u);
		this->forNeighborsOf(u, [&](node v) {
			if (marked[v] == 0) {
				q.push(v);
				marked[v] = 1;
			}
		});
	} while (!q.empty());
}

template<typename L>
inline void NetworKit::Graph::forEdgesOf(node u, L handle) {
	for (node v : this->adja[u]) {
		if (v != none) {
			handle(u, v);
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forEdgesOf(node u, L handle) const {
	for (node v : this->adja[u]) {
		if (v != none) {
			handle(u, v);
		}
	}
}

template<typename L>
void NetworKit::Graph::forEdgesOfInDegreeIncreasingOrder(node u, L handle) const {
	// TODO: iterating over neighbors ordered by degree does not need privileged access to graphs data structure and should
	// therefore be implemented inside the algorithm. - cls
	auto hasSmallerDegree = [&](node v1, node v2) {
		return degree(v1) < degree(v2); // FIXME
	};

	StdVector<node> neighbors = adja[u];
	std::sort(neighbors.begin(), neighbors.end(), hasSmallerDegree);

	for (node v : neighbors) {
		if (v != none) {
			handle(u, v);
		}
	}
}

template<typename L>
void NetworKit::Graph::forEdgesOfInDegreeIncreasingOrder(node u, L handle) {
	auto hasSmallerDegree = [&](node v1, node v2) {
		return degree(v1) < degree(v2); // FIXME
	};

	StdVector<node> neighbors = adja[u];
	std::sort(neighbors.begin(), neighbors.end(), hasSmallerDegree);

	for (node v : neighbors) {
		if (v != none) {
			handle(u, v);
		}
	}
}


template<typename L>
inline void NetworKit::Graph::parallelForNodePairs(L handle) {
#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				// call node pair function
				if (exists[v]) {
					handle(u, v);
				}
			}
		}

	}
}

template<typename L>
inline void NetworKit::Graph::parallelForNodePairs(L handle) const {
#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				// call node pair function
				if (exists[v]) {
					handle(u, v);
				}
			}
		}

	}
}

template<typename L>
inline void NetworKit::Graph::breadthFirstEdgesFrom(node r, L handle) {
	// TODO: implement BFS iterator for edges
	throw std::runtime_error("TODO");
}

template<typename L>
inline void NetworKit::Graph::forWeightedEdges(L handle) {
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				if (weighted) {
					edgeweight w = this->eweights[u][vi];
					handle(u, v, w);
				} else {
					handle(u, v, defaultEdgeWeight);
				}
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forWeightedEdges(L handle) const {
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				if (weighted) {
					edgeweight w = this->eweights[u][vi];
					handle(u, v, w);
				} else {
					handle(u, v, defaultEdgeWeight);
				}
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::parallelForWeightedEdges(L handle) {
#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				if (weighted) {
					edgeweight w = this->eweights[u][vi];
					handle(u, v, w);
				} else {
					handle(u, v, defaultEdgeWeight);
				}
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::parallelForWeightedEdges(L handle) const {
#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				if (weighted) {
					edgeweight w = this->eweights[u][vi];
					handle(u, v, w);
				} else {
					handle(u, v, defaultEdgeWeight);
				}
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forWeightedEdgesOf(node u, L handle) {
	const count asize = (count) adja[u].size();
	for (index i = 0; i < asize; ++i) {
		node v = adja[u][i];
		if (v != none) {
			if (weighted) {
				edgeweight w = this->eweights[u][i];
				handle(u, v, w);
			} else {
				handle(u, v, defaultEdgeWeight);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forWeightedEdgesOf(node u, L handle) const {
	for (index i = 0; i < adja[u].size(); ++i) {
		node v = adja[u][i];
		if (v != none) {
			if (weighted) {
				edgeweight w = this->eweights[u][i];
				handle(u, v, w);
			} else {
				handle(u, v, defaultEdgeWeight);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forNodesWithAttribute(std::string attrKey,
		L handle) {
	// get nodemap for attrKey

//	auto nodeMap; // ?
//
//	auto findIdPair = this->attrKey2IdPair.find(attrKey);
//	if (findIdPair != this->attrKey2IdPair.end()) {
//		std::pair<index, index> idPair = findIdPair->second;
//		index typeId = idPair.first;
//		index mapId = idPair.second;
//
//		// nodemaps are in a vector, one for each node attribute type int, double, NodeAttribute
//		switch (typeId) {
//		case 0:
//			nodeMap = this->nodeMapsInt[mapId];
//			break;
//		case 1:
//			nodeMap = this->nodeMapsdouble[mapId];
//			break;
//		}
//
//		// iterate over nodes and call handler with attribute
//		this->forNodes([&](node u) {
//			auto attr = nodeMap[u];
//			handle(u, attr);
//		});
//	} else {
//		throw std::runtime_error("node attribute not found");
//	}

	// TODO: iterate over nodes with atributes
	throw std::runtime_error("TODO");
}

template<typename L>
inline void NetworKit::Graph::forEdgesWithAttribute_double(int attrId,
		L handle) {
	std::vector<std::vector<double> > edgeMap = this->edgeMaps_double[attrId];
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < (index) adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			double attr = edgeMap[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				handle(u, v, attr);
			}
		}
	}
}

template<typename L>
inline void NetworKit::Graph::forEdgesWithAttribute_double(int attrId,
		L handle) const {
	std::vector<std::vector<double> > edgeMap = this->edgeMaps_double[attrId];
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < (index) adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			double attr = edgeMap[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				handle(u, v, attr);
			}
		}
	}
}

template<typename C, typename L>
inline void NetworKit::Graph::forNodesWhile(C condition, L handle) {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			if (!condition()) {
				break; // if condition does not hold, break from loop and do not call handle
			}
			handle(v);

		}
	}
}



template<typename C, typename L>
inline void NetworKit::Graph::forNodes(C condition, L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			if (!condition()) {
				break; // if condition does not hold, break from loop and do not call handle
			}
			handle(v);
		}
	}
}

template<typename L>
void NetworKit::Graph::forNodesInRandomOrder(L handle) {
	std::vector<node> randVec(z);
	for (node v = 0; v < z; ++v) {
		randVec[v] = v;
	}
	random_shuffle(randVec.begin(), randVec.end());

	for (node v = 0; v < z; ++v) {
		node randv = randVec[v];
		if (exists[randv]) {
			handle(randv);
		}
	}
}

template<typename L>
void NetworKit::Graph::forNodesInRandomOrder(L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}


template<typename L>
void NetworKit::Graph::BFSfrom(node r, L handle) {
	std::vector<bool> marked(z);
	std::queue<node> q;
	q.push(r); // enqueue root
	marked[r] = true;
	do {
		node u = q.front();
		q.pop();
		// apply function
		handle(u);
		this->forNeighborsOf(u, [&](node v) {
			if (!marked[v]) {
				q.push(v);
				marked[v] = true;
			}
		});
	} while (!q.empty());
};


template<typename L>
void NetworKit::Graph::DFSfrom(node r, L handle) {
	std::vector<bool> marked(z);
	std::stack<node> stack;
	stack.push(r); // enqueue root
	marked[r] = true;
	do {
		node u = stack.top();
		stack.pop();
		// apply function
		handle(u);
		this->forNeighborsOf(u, [&](node v) {
			if (!marked[v]) {
				stack.push(v);
				marked[v] = true;
			}
		});
	} while (!stack.empty());
};



#endif /* GRAPH_H_ */
