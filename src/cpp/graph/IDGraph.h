/*
 * IDGraph.h
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef IDGRAPH_H
#define IDGRAPH_H

#include "IGraph.h"

namespace NetworKit {

/**
 * Interface for all directed graph classes. Every graph class has to implement all interface methods.
 */
class IDGraph : public virtual IGraph {

public:

	IDGraph() = default;

	/* move assignments and more currently removed because a problems with gcc 4.7.1 on phipute1 */
	// IDGraph(const IDGraph& other) = default;
	// IDGraph(IDGraph&& other) = default;
	// IDGraph& operator=(IDGraph&& other) = default;
	// IDGraph& operator=(const IDGraph& other) = default;

	/** NODE PROPERTIES **/

	/**
	 * @return true if the node is isolated (= degree is 0)
	 */
	bool isIsolated(node v) const { return degreeIn(v) == 0 && degreeOut(v) == 0; }

	/**
	 * Return the number of neighbors for node v. For a Directed Graph return the number of outgoing neighbors for node v.
	 */
	virtual count degree(node v) const { return degreeOut(v); }

	/**
	 * Return the number of incoming edges to node v.
	 */
	virtual count degreeIn(node v) const = 0;

	/**
	 * Return the number of outgoing edges from node v.
	 */
	virtual count degreeOut(node v) const = 0;

	/**
	 * @return Weighted degree of @a v. For directed graphs this is the sum of weights off all outgoing edges fo @a v.
	 */
	edgeweight weightedDegree(node v) const { return weightedDegreeOut(v); }

	/**
	 * @return Weighted degree of @a v.
	 */
	virtual edgeweight weightedDegreeIn(node v) const = 0;

	/**
	 * @return Weighted degree of @a v.
	 */
	virtual edgeweight weightedDegreeOut(node v) const = 0;


	/** NEIGHBORHOOD ITERATORS **/

	/**
	 * Iterate over all neighbors of a node and call handler (lamdba closure).
	 */
	virtual void forNeighborsOf(node u, FNode f) const { forOutNeighborsOf(u, f); }

	/**
	 * Iterate over all adjacent nodes, which have an edge from u.
	 */
	virtual void forOutNeighborsOf(node u, FNode f) const = 0;

	/**
	 * Iterate over all adjacent nodes, which have an edge to u.
	 */
	virtual void forInNeighborsOf(node u, FNode f) const = 0;

	/**
	 * Iterate over all edge weights of a node and call handler (lamdba closure).
	 */
	virtual void forWeightedNeighborsOf(node u, FNodeWeighted f) const { forWeightedOutNeighborsOf(u, f); }

	/**
	 * Iterate over all outgoing edge weights of a node and call handler (lamdba closure).
	 */
	virtual void forWeightedOutNeighborsOf(node u, FNodeWeighted f) const = 0;
	
	/**
	 * Iterate over all incoming edge weights of a node and call handler (lamdba closure).
	 */
	virtual void forWeightedInNeighborsOf(node u, FNodeWeighted f) const = 0;

	/**
	 * Iterate over all incident edges of a node and call handler (lamdba closure).
	 */
	virtual void forEdgesOf(node u, FEdge f) const { forOutEdgesOf(u, f); }

	/**
	 * Iterate over all outgoing edges of the graph and call handler (lambda closure).
	 */
	virtual void forOutEdgesOf(node u, FEdge f) const = 0;
	
	/**
	 * Iterate over all incoming edges of the graph and call handler (lambda closure).
	 */
	virtual void forInEdgesOf(node u, FEdge f) const = 0;

	/**
	 * Iterate over all incident edges of a node and call handler (lamdba closure).
	 *
	 * Handle takes parameters (u, v, w) where w is the edge weight.
	 *
	 */
	virtual void forWeightedEdgesOf(node u, FEdgeWeighted f) const { forWeightedOutEdgesOf(u, f); }
	
	/**
	 * Iterate over all incident edges from u and call handler (lamdba closure).
	 *
	 * Handle takes parameters (u, v, w) where w is the edge weight.
	 *
	 */
	virtual void forWeightedOutEdgesOf(node u, FEdgeWeighted f) const = 0;

	/**
	 * Iterate over all incident edges from u and call handler (lamdba closure).
	 *
	 * Handle takes parameters (u, v, w) where w is the edge weight.
	 *
	 */
	virtual void forWeightedInEdgesOf(node u, FEdgeWeighted f) const = 0;

};

} /* namespace NetworKit */

#endif /* IDGRAPH_H */
