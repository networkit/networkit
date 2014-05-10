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

#include "../Globals.h"
#include "../viz/Point.h"

namespace NetworKit {

/** Typedefs **/

typedef uint64_t index; // more expressive name for an index into an array
typedef uint64_t count; // more expressive name for an integer quantity
typedef index node; // node indices are 0-based
typedef double edgeweight; // edge weight type

constexpr index none = std::numeric_limits<index>::max();

/**
 * Interface for all graph classes. Every graph class has to implement all interface methods.
 */
class IGraph {

protected:
	count graphId;

public:
	IGraph() {
		static uint64_t nextGraphId = 1;
		graphId = nextGraphId++;
	}

	/**
	 * Get the ID of this graph. The ID is a unique unsigned integer given to
	 * every graph on construction.
	 */
	uint64_t getId() const { return graphId; }

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
	 * Return the number of neighbors for node v. For directed graphs this is the sum of
	 * in- and outgoing edges.
	 */
	virtual count degree(node v) const = 0;

	/**
	 * @return random node of the graph
	 */
	virtual node randomNode() const = 0;


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


	/** NODE ITERATORS **/

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodes(L handle) {};

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodes(L handle) const {};

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodes(L handle) {};

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodes(L handle) const {};

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure) as long as the condition remains true.
	 * This allows for breaking from a node loop.
	 */
	template<typename C, typename L> void forNodesWhile(C condition, L handle) {};

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure) as long as the condition remains true.
	 * This allows for breaking from a node loop.
	 */
	template<typename C, typename L> void forNodes(C condition, L handle) const {};

	/**
	 * Iterate randomly over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodesInRandomOrder(L handle) {};

	/**
	 * Iterate randomly over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodesInRandomOrder(L handle) const {};

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 * Using schedule(guided) to remedy load-imbalances due to e.g. unequal degree distribution.
	 */
	template<typename L> void balancedParallelForNodes(L handle) {};

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 * Using schedule(guided) to remedy load-imbalances due to e.g. unequal degree distribution.
	 */
	template<typename L> void balancedParallelForNodes(L handle) const {};

	/**
	 * Iterate over all undirected pairs of nodesand call handler (lambda closure).
	 */
	template<typename L> void forNodePairs(L handle) {};

	/**
	 * Iterate over all undirected pairs of nodesand call handler (lambda closure).
	 */
	template<typename L> void forNodePairs(L handle) const {};

	/**
	 * Iterate over all undirected pairs of nodes in parallel and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodePairs(L handle) {};

	/**
	 * Iterate over all undirected pairs of nodes in parallel and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodePairs(L handle) const {};


 	/** REDUCTION ITERATORS **/

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForNodes(L handle) { return 0.0; };

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForNodes(L handle) const { return 0.0; };


	/** EDGE ITERATORS **/

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void forEdges(L handle) {};

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void forEdges(L handle) const {};

	/**
	 * Iterate in parallel over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForEdges(L handle) {};

	/**
	 * Iterate in parallel over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForEdges(L handle) const {};

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 *
	 * Handler takes arguments (u, v, w) where u and v are the nodes of the edge and w is its weight.
	 *
	 */
	template<typename L> void forWeightedEdges(L handle) {};

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
	 *
	 */
	template<typename L> void parallelForWeightedEdges(L handle) {};

	/**
	 * Iterate over all edges of the graph and call handler (lambda closure).
	 *
	 * Handler takes arguments (u, v, w) where u and v are the nodes of the edge and w is its weight.
	 */
	template<typename L> void parallelForWeightedEdges(L handle) const {};
};





} /* namespace NetworKit */

#endif /* IGRAPH_H */
