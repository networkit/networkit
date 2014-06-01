/*
 * AbstractGraph.h
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef ABSTRACT_GRAPH_H
#define ABSTRACT_GRAPH_H

#include <vector>
#include <queue>
#include <stack>
// #include <tbb/concurrent_vector.h>

#include "IGraph.h"
#include "Coordinates.h"
#include "../Globals.h"

namespace NetworKit {

// TODO: test tbb::concurrent_vector
// template <typename T> using _Vector = std::vector<T>;

/**
 * An abstract graph class providing general implementations of directed and undirected graphs (with optional weights).
 */
class AbstractGraph : public virtual IGraph {

protected:

	// scalars
	count n; //!< current number of nodes
	count m; //!< current number of edges
	node z; //!< current upper bound of node ids
	count t; //!< current time step
	bool weighted; //!< true if this graph supports edge weights other than 1.0

	// per node data
	// Vector<count> deg; //!< degree of each node (size of neighborhood)
	std::vector<bool> exists; //!< exists[v] is true if node v has not been removed from the graph
	Coordinates<float> coordinates; //!< coordinates of nodes (if present)

	// per edge data
	// Vector<Vector<node> > adja; //!< neighbors/adjacencies
	// Vector<Vector<edgeweight> > eweights; //!< edge weights

	// graph attributes
	std::string name;

	// user-defined edge attributes

	//	attribute maps storage

	// std::vector<std::vector<std::vector<double> > > edgeMaps_double; // contains edge maps (u, v) -> double

	// defaults

	// std::vector<double> edgeAttrDefaults_double; // stores default value for edgeMaps_double[i] at index i


public:

	// defaults
	static constexpr double defaultEdgeWeight = 1.0;
	static constexpr edgeweight nullWeight = 0.0;

	/** ATTRIBUTE ABSTRACT BASE CLASSES **/

	// class NodeAttribute {
	// 	// abstract
	// };

	// class EdgeAttribute {
	// 	// abstract
	// };


	/** GRAPH INTERFACE **/

	/** 
	 * Create a graph of n nodes.
	 */
	AbstractGraph(count n = 0, bool weighted = false);

	/* move assignments and more currently removed because a problems with gcc 4.7.1 on phipute1 */
	AbstractGraph(const AbstractGraph& other) = default;
	AbstractGraph(AbstractGraph&& other) = default;
	AbstractGraph& operator=(AbstractGraph&& other) = default;
	AbstractGraph& operator=(const AbstractGraph& other) = default;

	/**
	 * Calculate an approximation of the memory used by this graph. Only memory increasing with the
	 * number of edges or nodes of this graph is taken into account. 
	 */
	virtual count getMemoryUsage() const ;

	/**
	 * Try to save some memory by shrinking internal data structures of the graph. Only run this
	 * once you finished editing the graph. Otherwise it will cause unnecessary reallocation of
	 * memory. 
	 */
	virtual void shrinkToFit();

	/**
	 * Set name of graph.
	 */
	virtual void setName(std::string name);

	/*
	 * @return name of graph
	 */
	virtual std::string getName() const;

	/**
	 * Get string representation
	 */
	virtual std::string toString() const;


	/** NODE MODIFIERS **/

	/**
	 * Add a new node to the graph and return it.
	 */
	virtual node addNode();

	/**
	 * Add a new node to the graph with coordinates @a x and @y and return it.
	 */
	virtual node addNode(float x, float y);

	/**
	 * Remove an isolated node from the graph.
	 *
	 * Although it would be convenient to remove all incident edges at the same time,
	 * this causes complications for dynamic applications. Therefore, removeNode is an
	 * atomic event. All incident edges need to be removed first and an exception is thrown
	 * otherwise.
	 */
	virtual void removeNode(node u);

	/**
	 * Check if node exists in the graph.
	 */
	virtual bool hasNode(node u) const;


	/** NODE PROPERTIES **/

	/**
	 * @return random node of the graph
	 */
	virtual node randomNode() const;


	/** GLOBAL PROPERTIES **/

	/**
	 * Return true if this graph supports edge weights other than 1.0
	 */
	virtual bool isWeighted() const;

	/** 
	 * Return true if this graph supports directed edges.
	 */
	virtual bool isDirected() const = 0;

	/**
	 * Return true if graph contains no nodes.
	 */
	virtual bool isEmpty() const;

	/**
	 * Return the number of nodes in the graph.
	 *
	 */
	virtual count numberOfNodes() const;

	/**
	 * Return the number of edges in the graph.
	 *
	 *	 */
	virtual count numberOfEdges() const;

	/**
	 * @return the number of loops {v, v} in the graph.
	 *
	 * This involves calculation, so store result if needed multiple times.
	 */
	virtual count numberOfSelfLoops() const;

 	/**
	 * Get an upper bound for the node ids in the graph.
	 */
	virtual index upperNodeIdBound() const;


	/** DYNAMICS **/

	/**
	 * Trigger a time step - increments counter.
	 */
	virtual void timeStep();

	/**
	 * Get time step counter.
	 */
	virtual count time();


	/** COORDINATES **/

	virtual void setCoordinate(node v, Point<float> value);

	virtual Point<float>& getCoordinate(node v);

	virtual float minCoordinate(count dim);

	virtual float maxCoordinate(count dim);

	virtual void initCoordinates();

	/** SUMS **/
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
	 * Iterate over all nodes of the graph and call fr (lambda closure).
	 */
	virtual void forNodes(FNode f) const;

	/**
	 * Iterate in parallel over all nodes of the graph and call fr (lambda closure).
	 */
	virtual void parallelForNodes(FNode f) const;

	/**
	 * Iterate over all nodes of the graph and call fr (lambda closure) as long as the condition remains true.
	 * This allows for breaking from a node loop.
	 */
	virtual void forNodesWhile(FCondition condition, FNode f) const;

	/**
	 * Iterate over all nodes of the graph and call fr (lambda closure) as long as the condition remains true.
	 * This allows for breaking from a node loop.
	 */
	virtual void forNodes(FNodeCondition condition, FNode f) const;

	/**
	 * Iterate randomly over all nodes of the graph and call fr (lambda closure).
	 */
	virtual void forNodesInRandomOrder(FNode f) const;

	/**
	 * Iterate in parallel over all nodes of the graph and call fr (lambda closure).
	 * Using schedule(guided) to remedy load-imbalances due to e.g. unequal degree distribution.
	 */
	virtual void balancedParallelForNodes(FNode f) const;

	/**
	 * Iterate over all undirected pairs of nodesand call fr (lambda closure).
	 */
	virtual void forNodePairs(FNodePair f) const;

	/**
	 * Iterate over all undirected pairs of nodes in parallel and call fr (lambda closure).
	 */
	virtual void parallelForNodePairs(FNodePair f) const;


 	/** REDUCTION ITERATORS **/

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the fr
	 */
	virtual double parallelSumForNodes(FNodeSum f) const;


	/** GRAPH SEARCHES **/

	virtual void BFSfrom(node r, FNode f) const;

	virtual void BFSEdgesfrom(node r, FEdge f) const;

	virtual void DFSfrom(node r, FNode f) const;
	
	virtual void DFSEdgesfrom(node r, FEdge f) const;

};

} /* namespace NetworKit */


/** NODE ITERATORS **/

inline void NetworKit::AbstractGraph::forNodes(NetworKit::FNode f) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			f(v);
		}
	}
}

inline void NetworKit::AbstractGraph::parallelForNodes(NetworKit::FNode f) const {
	#pragma omp parallel for
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			f(v);
		}
	}
}

inline void NetworKit::AbstractGraph::forNodesWhile(NetworKit::FCondition condition, NetworKit::FNode f) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			if (!condition()) {
				break; // if condition does not hold, break from loop and do not call f
			}
			f(v);

		}
	}
}

inline void NetworKit::AbstractGraph::forNodes(NetworKit::FNodeCondition condition, NetworKit::FNode f) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			if (!condition(v)) {
				f(v);
			}
		}
	}
}

inline void NetworKit::AbstractGraph::forNodesInRandomOrder(NetworKit::FNode f) const {
	std::vector<node> randVec = nodes();
	random_shuffle(randVec.begin(), randVec.end());
	for (node v : randVec) {
		f(v);
	}
}

inline void NetworKit::AbstractGraph::balancedParallelForNodes(NetworKit::FNode f) const {
	#pragma omp parallel for schedule(guided) // TODO: define min block size (and test it!)
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			f(v);
		}
	}
}

inline void NetworKit::AbstractGraph::forNodePairs(NetworKit::FNodePair f) const {
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				if (exists[v]) {
					f(u, v);
				}
			}
		}

	}
}

inline void NetworKit::AbstractGraph::parallelForNodePairs(NetworKit::FNodePair f) const {
	#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				if (exists[v]) {
					f(u, v);
				}
			}
		}
	}
}


/** REDUCTION ITERATORS **/

inline double NetworKit::AbstractGraph::parallelSumForNodes(NetworKit::FNodeSum f) const {
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			sum += f(v);
		}
	}
	return sum;
}

#endif /* ABSTRACT_GRAPH_H */
