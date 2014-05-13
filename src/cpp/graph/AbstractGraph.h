/*
 * AbstractGraph.h
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef ABSTRACT_GRAPH_H
#define ABSTRACT_GRAPH_H

// #include <functional>
// #include <cassert>
#include <vector>
// #include <cinttypes>
// #include <queue>
// #include <stack>
// #include <map>
// #include <set>
// #include <sstream>
// #include <cstdint>
// #include <algorithm>
// #include <tbb/concurrent_vector.h>

#include "IGraph.h"
#include "Coordinates.h"
// #include "../auxiliary/Log.h"
#include "../Globals.h"


namespace NetworKit {

// TODO: test tbb::concurrent_vector
template <typename T> using _Vector = std::vector<T>;


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
	_Vector<bool> exists; //!< exists[v] is true if node v has not been removed from the graph
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

	AbstractGraph(const AbstractGraph& other) = default;

	AbstractGraph(AbstractGraph&& other) = default;

	~AbstractGraph() = default;

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
	 * Iterate over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodes(L handle);

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodes(L handle) const;

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodes(L handle);

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodes(L handle) const;

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure) as long as the condition remains true.
	 * This allows for breaking from a node loop.
	 */
	template<typename C, typename L> void forNodesWhile(C condition, L handle);

	/**
	 * Iterate over all nodes of the graph and call handler (lambda closure) as long as the condition remains true.
	 * This allows for breaking from a node loop.
	 */
	template<typename C, typename L> void forNodes(C condition, L handle) const;

	/**
	 * Iterate randomly over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodesInRandomOrder(L handle);

	/**
	 * Iterate randomly over all nodes of the graph and call handler (lambda closure).
	 */
	template<typename L> void forNodesInRandomOrder(L handle) const;

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 * Using schedule(guided) to remedy load-imbalances due to e.g. unequal degree distribution.
	 */
	template<typename L> void balancedParallelForNodes(L handle);

	/**
	 * Iterate in parallel over all nodes of the graph and call handler (lambda closure).
	 * Using schedule(guided) to remedy load-imbalances due to e.g. unequal degree distribution.
	 */
	template<typename L> void balancedParallelForNodes(L handle) const;

	/**
	 * Iterate over all undirected pairs of nodesand call handler (lambda closure).
	 */
	template<typename L> void forNodePairs(L handle);

	/**
	 * Iterate over all undirected pairs of nodesand call handler (lambda closure).
	 */
	template<typename L> void forNodePairs(L handle) const;

	/**
	 * Iterate over all undirected pairs of nodes in parallel and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodePairs(L handle);

	/**
	 * Iterate over all undirected pairs of nodes in parallel and call handler (lambda closure).
	 */
	template<typename L> void parallelForNodePairs(L handle) const;


 	/** REDUCTION ITERATORS **/

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForNodes(L handle);

	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForNodes(L handle) const;

};

} /* namespace NetworKit */

/** NODE ITERATORS **/

template<typename L>
inline void NetworKit::AbstractGraph::forNodes(L handle) {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::AbstractGraph::forNodes(L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::AbstractGraph::parallelForNodes(L handle) {
	#pragma omp parallel for
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::AbstractGraph::parallelForNodes(L handle) const {
	#pragma omp parallel for
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename C, typename L>
inline void NetworKit::AbstractGraph::forNodesWhile(C condition, L handle) {
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
inline void NetworKit::AbstractGraph::forNodes(C condition, L handle) const {
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
void NetworKit::AbstractGraph::forNodesInRandomOrder(L handle) {
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
void NetworKit::AbstractGraph::forNodesInRandomOrder(L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::AbstractGraph::balancedParallelForNodes(L handle) {
	#pragma omp parallel for schedule(guided) // TODO: define min block size (and test it!)
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::AbstractGraph::balancedParallelForNodes(L handle) const {
	#pragma omp parallel for schedule(guided) // TODO: define min block size (and test it!)
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<typename L>
inline void NetworKit::AbstractGraph::forNodePairs(L handle) {
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
inline void NetworKit::AbstractGraph::forNodePairs(L handle) const {
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
inline void NetworKit::AbstractGraph::parallelForNodePairs(L handle) {
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

template<typename L>
inline void NetworKit::AbstractGraph::parallelForNodePairs(L handle) const {
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


/** REDUCTION ITERATORS **/

template<typename L>
inline double NetworKit::AbstractGraph::parallelSumForNodes(L handle) {
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
inline double NetworKit::AbstractGraph::parallelSumForNodes(L handle) const {
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			sum += handle(v);
		}
	}
	return sum;
}


#endif /* ABSTRACT_GRAPH_H */
