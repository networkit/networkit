/*
 * IGraph.h
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef IGRAPH_H
#define IGRAPH_H

#include <limits>
#include <string>
#include <algorithm>
#include <stdexcept>

#include "../Globals.h"
#include "../viz/Point.h"

#include "BasicGraph.h"

namespace NetworKit {

/** Typedefs **/
// moved to Globals.h

// typedef uint64_t index; // more expressive name for an index into an array
// typedef uint64_t count; // more expressive name for an integer quantity
// typedef index node; // node indices are 0-based
// typedef double edgeweight; // edge weight type

// constexpr index none = std::numeric_limits<index>::max();

/**
 * Interface for all graph classes. Every graph class has to implement all interface methods.
 */
class IGraph {

protected:
	count graphId;

public:

	IGraph() {
		static count nextGraphId = 1;
		graphId = nextGraphId++;
	}

	/* move assignments and more currently removed because a problems with gcc 4.7.1 on phipute1 */
	// IGraph(const IGraph& other) = default;
	// IGraph(IGraph&& other) = default;
	// IGraph& operator=(IGraph&& other) = default;
	// IGraph& operator=(const IGraph& other) = default;

	/**
	 * Get the ID of this graph. The ID is a unique unsigned integer given to
	 * every graph on construction.
	 */
	count getId() const { return graphId; }

	/**
	 * Calculate an approximation of the memory used by this graph. Only memory increasing with the
	 * number of edges or nodes of this graph is taken into account. 
	 */
	virtual count getMemoryUsage() const = 0;

	/**
	 * Try to save some memory by shrinking internal data structures of the graph. Only run this
	 * once you finished editing the graph. Otherwise it will cause unnecessary reallocation of
	 * memory. 
	 */
	virtual void shrinkToFit() = 0;

	/**
	 * Set name of graph.
	 */
	virtual void setName(std::string name) = 0;

	/*
	 * @return name of graph
	 */
	virtual std::string getName() const = 0;

	/**
	 * Get string representation
	 */
	virtual std::string toString() const = 0;


	/** NODE MODIFIERS **/

	/**
	 * Add a new node to the graph and return it.
	 */
	virtual node addNode() = 0;

	/**
	 * Add a new node to the graph with coordinates @a x and @y and return it.
	 */
	virtual node addNode(float x, float y) = 0;

	/**
	 * Remove an isolated node v from the graph.
	 *
	 * Although it would be convenient to remove all incident edges at the same time,
	 * this causes complications for dynamic applications. Therefore, removeNode is an
	 * atomic event. All incident edges need to be removed first and an exception is thrown
	 * otherwise.
	 */
	virtual void removeNode(node v) = 0;

	/**
	 * Check if node v exists in the graph.
	 */
	virtual bool hasNode(node v) const = 0;


	/** NODE PROPERTIES **/

	/**
	 * Return the number of neighbors for node v.
	 */
	virtual count degree(node v) const = 0;

	/**
	 * @return true if the node is isolated (= degree is 0)
	 */
	virtual bool isIsolated(node v) const = 0;

	/**
	 * @return Weighted degree of @a v. For directed graphs this is the sum of weights off all outgoing edges fo @a v.
	 */
	virtual edgeweight weightedDegree(node v) const = 0;

	/**
	 * @return random node of the graph
	 */
	virtual node randomNode() const = 0;


	/** EDGE MODIFIERS **/

	/**
	 * Insert an directed edge between from @a u to @a v.
	 */
	virtual void addEdge(node u, node v, edgeweight weight) = 0;

	/**
	 * Remove directed edge between from @a u to @a v.
	 */
	virtual void removeEdge(node u, node v) = 0;

	/**
	 * Check if directed edge {u,v} exists.
	 *
	 */
	virtual bool hasEdge(node u, node v) const = 0;

	/**
	 * Merges edge {u,v} to become a supernode. Edges to u and v are
	 * rewired, multiple edges merged and their weights added.
	 * The vertex weights of @a u and @a v are added.
	 * A self-loop is only created if @a discardSelfLoop is set to false.
	 *
	 * @return New node that has been created if u != v. Otherwise none.
	 */
	// node mergeEdge(node u, node v, bool discardSelfLoop = true);


	/** GLOBAL PROPERTIES **/

	/**
	 * Return true if this graph supports edge weights other than 1.0
	 */
	virtual bool isWeighted() const = 0;

	/** 
	 * Return true if this graph supports directed edges.
	 */
	virtual bool isDirected() const = 0;

	/**
	 * Return true if graph contains no nodes.
	 */
	virtual bool isEmpty() const = 0;

	/**
	 * Return the number of nodes in the graph.
	 *
	 */
	virtual count numberOfNodes() const = 0;

	/**
	 * Return the number of edges in the graph.
	 *
	 *	 */
	virtual count numberOfEdges() const = 0;

	/**
	 * @return the number of loops {v, v} in the graph.
	 *
	 * This involves calculation, so store result if needed multiple times.
	 */
	virtual count numberOfSelfLoops() const = 0;

 	/**
	 * Get an upper bound for the node ids in the graph.
	 */
	virtual index upperNodeIdBound() const = 0;

	/** DYNAMICS **/

	/**
	 * Trigger a time step - increments counter.
	 */
	virtual void timeStep() = 0;

	/**
	 * Get time step counter.
	 */
	virtual count time() = 0;


	/** COORDINATES **/

	virtual void setCoordinate(node v, Point<float> value) = 0;

	virtual Point<float>& getCoordinate(node v) = 0;

	virtual float minCoordinate(count dim) = 0;

	virtual float maxCoordinate(count dim) = 0;

	virtual void initCoordinates() = 0;


	/** EDGE ATTRIBUTES **/

	/**
	 * Return edge weight.
	 *
	 * Return 0 if edge does not exist.
	 */
	virtual edgeweight weight(node u, node v) const = 0;

	/**
	 * Set the weight of an edge. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	virtual void setWeight(node u, node v, edgeweight w) = 0;

	/**
	 * Increase the weight of an edge. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	virtual void increaseWeight(node u, node v, edgeweight w) = 0;

	/**
	 * Add new edge map for an attribute of type double.
	 */
	virtual int addEdgeAttribute_double(double defaultValue) = 0;

	/**
	 * @return attribute of type double for an edge.
	 *
	 * @param[in]	u	node
	 * @param[in]	v	node
	 * @param[in]	attrId	attribute id
	 */
	virtual double attribute_double(node u, node v, int attrId) const = 0;

	/**
	 * Set edge attribute of type double If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	attr	double edge attribute
	 */
	virtual void setAttribute_double(node u, node v, int attrId, double attr) = 0;


	/** SUMS **/

	/**
	 * @return sum of all edge weights
	 */
	virtual edgeweight totalEdgeWeight() const = 0;


	/** Collections **/

	/**
	 * Return list of nodes
	 */
	virtual std::vector<node> nodes() const = 0;

	/**
	 * Return list of edges as node pairs.
	 */
	virtual std::vector<std::pair<node, node> > edges() const = 0;


	/**
	 * Return list of neighbors for given node.
	 */
	virtual std::vector<node> neighbors(node u) const = 0;


	/** NODE ITERATORS **/

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure).
	 */
	// template<typename L> void forNodes(L handle) const {};
	virtual void forNodes(FNode f) const = 0;

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodes(L handle) const {};

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure) as long as the condition remains true.
	 * This allows for breaking from a node loop.
	 */
	template<typename C, typename L> void forNodesWhile(C condition, L handle) const {};

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure) as long as the condition remains true.
	 * This allows for breaking from a node loop.
	 */
	template<typename C, typename L> void forNodes(C condition, L handle) const {};

	/**
	 * Iterate randomly over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodesInRandomOrder(L handle) const {};

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 * Using schedule(guided) to remedy load-imbalances due to e.g. unequal degree distribution.
	 */
	template<typename L> void balancedParallelForNodes(L handle) const {};

	/**
	 * Iterate over all undirected pairs of nodesand call handler (lambda closure).
	 */
	template<typename L> void forNodePairs(L handle) const {};

	/**
	 * Iterate over all undirected pairs of nodes in parallel and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodePairs(L handle) const {};


 	/** REDUCTION ITERATORS **/

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForNodes(L handle) const { return 0.0; };


	/** EDGE ITERATORS **/

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 */
	// template<typename L> void forEdges(L handle) const {};
	virtual void forEdges(FEdge f) const = 0;

	/**
	 * Iterate in parallel over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForEdges(L handle) const {};

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 *
	 * Handler takes arguments (u, v, w) where u and v are the nodes of the edge and w is its weight.
	 */
	template<typename L> void forWeightedEdges(L handle) const {};

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 *
	 * Handler takes arguments (u, v, w) where u and v are the nodes of the edge and w is its weight.
	 */
	template<typename L> void parallelForWeightedEdges(L handle) const {};

	/**
	 * Iterate over all edges of the const graph and call handler (lambda closure).
	 *
	 *	@param[in]	attrId		attribute id
	 *	@param[in]	handle 		takes arguments (u, v, a) where a is an edge attribute of edge {u, v}
	 *
	 */
	template<typename L> void forEdgesWithAttribute_double(int attrId, L handle) const {};


	/** NEIGHBORHOOD ITERATORS **/

	/**
	 * Iterate over all neighbors of a node and call handler (lamdba closure).
	 */
	// template<typename L> void forNeighborsOf(node u, L handle) const {}
	virtual void forNeighborsOf(node u, FNode handle) const = 0;

	/**
	 * Iterate over all edge weights of a node and call handler (lamdba closure).
	 */
	template<typename L> void forWeightedNeighborsOf(node u, L handle) const {}

	/**
	 * Iterate over all incident edges of a node and call handler (lamdba closure).
	 */
	template<typename L> void forEdgesOf(node u, L handle) const {}

	/**
	 * Iterate over all incident edges of a node and call handler (lamdba closure).
	 *
	 * Handle takes parameters (u, v, w) where w is the edge weight.
	 *
	 */
	template<typename L> void forWeightedEdgesOf(node u, L handle) const {}


	/** REDUCTION ITERATORS **/

	/**
	 * Iterate in parallel over all edges and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForWeightedEdges(L handle) const;


	/** GRAPH SEARCHES **/

	virtual void BFSfrom(node r, FNode f) const = 0;

	virtual void BFSEdgesfrom(node r, FEdge f) const = 0;

	virtual void DFSfrom(node r, FNode f) const = 0;
	
	virtual void DFSEdgesfrom(node r, FEdge f) const = 0;

	// template<typename L> void BFSfrom(node r, L handle) const;

	// template<typename L> void BFSEdgesfrom(node r, L handle) const;

	// template<typename L> void DFSfrom(node r, L handle) const;

	// template<typename L> void DFSEdgesfrom(node r, L handle) const;

};

} /* namespace NetworKit */

#endif /* IGRAPH_H */
