/*
 * BasicGraph.h
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef BASICGRAPH_H_
#define BASICGRAPH_H_

#include <algorithm>
#include <functional>
#include <vector>
#include <type_traits>
#include <utility>

#include "../Globals.h"
#include "Coordinates.h"
#include "../viz/Point.h"
#include "GraphData.h"

namespace NetworKit {

// typedef std::function<void(node)> FNode;
// typedef std::function<void(node, node)> FNodePair;
// typedef std::function<void(node, node)> FEdge;
// typedef std::function<void(node, node, double)> FEdgeWithWeight;
// typedef std::function<bool()> FCondition;

enum class Weighted {
	weighted,
	unweighted
};

enum class Directed {
	directed,
	undirected
};

// hide implementation details in extra namespace
namespace graph_impl {

// choose between types
template<bool b, class T1, class T2>
using either = typename std::conditional<b, T1, T2>::type;

// class for all graphs in NetworKit, some methods have special implementation depending on the template parameters
template<Weighted weighted, Directed directed>
class BasicGraph :
	private either<weighted == Weighted::weighted, WeightedData, UnweightedData>,
	private either<directed == Directed::directed, DirectedData, UndirectedData>
{
private:

	using DData = either<directed == Directed::directed, DirectedData, UndirectedData>;
	using WData = either<weighted == Weighted::weighted, WeightedData, UnweightedData>;

	// graph attributes
	count graphId;
	std::string name;

	// scalars
	count n; //!< current number of nodes
	count m; //!< current number of edges
	node z; //!< current upper bound of node ids
	count t; //!< current time step

	// per node data
	std::vector<bool> exists; //!< exists[v] is true if node v has not been removed from the graph
	Coordinates<float> coordinates; //!< coordinates of nodes (if present)
	std::vector<std::vector<edgeweight> > edgeWeights;
	// user-defined edge attributes
	// attribute maps storage
	std::vector<std::vector<std::vector<double> > > edgeMaps_double; // contains edge maps (u, v) -> double
	// default values
	std::vector<double> edgeAttrDefaults_double; // stores default value for edgeMaps_double[i] at index i

	index find(node u, node v) const { 
		throw std::runtime_error("find is deprecated, use indexInEdgeArray instead");
	}

public:
	BasicGraph(count n = 0);

	BasicGraph(const BasicGraph<weighted, directed>& other) = default;

	BasicGraph(BasicGraph<weighted, directed>&& other) = default;

	~BasicGraph() = default;

	BasicGraph<weighted, directed>& operator=(BasicGraph<weighted, directed>&& other) = default;

	BasicGraph<weighted, directed>& operator=(const BasicGraph<weighted, directed>& other) = default;


	/** GRAPH INFORMATION **/

	/**
	 * Get the ID of this graph. The ID is a unique unsigned integer given to
	 * every graph on construction.
	 */
	count getId() const { return graphId; }

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
	count getMemoryUsage() const;

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
	count degree(node v) const { return degree_impl(*this, v); }

	template<Weighted w>
	friend count degree_impl(const BasicGraph<w, directed>& G, node v);

	/**
	 * @return true if the node is isolated (= degree is 0)
	 */
	bool isIsolated(node v) const { return isIsolated(*this, v); }

	template<Weighted w>
	friend bool isIsolated_impl(const BasicGraph<w, directed>& G, node v);

	/**
	 * @return Weighted degree of @a v. For directed graphs this is the sum of weights off all outgoing edges fo @a v.
	 */
	edgeweight weightedDegree(node v) const;

	/**
	 * @return random node of the graph
	 */
	node randomNode() const;


	/** EDGE MODIFIERS **/

	/**
	 * Insert an directed edge between from @a u to @a v.
	 */
	void addEdge(node u, node v, edgeweight ew = defaultEdgeWeight) { addEdge_impl(*this, u, v, ew); }

	template<Weighted w>
	friend void addEdge_impl(BasicGraph<w, directed>& G, node u, node v, edgeweight ew);

	/**
	 * Remove directed edge between from @a u to @a v.
	 */
	void removeEdge(node u, node v) { removeEdge_impl(*this, u, v); }

	template<Weighted w>
	friend void removeEdge_impl(BasicGraph<w, directed>& G, node u, node v);

	/**
	 * Check if directed edge {u,v} exists.
	 *
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

	node mergeEdge(node u, node v, bool discardSelfLoop = true);


	/** GLOBAL PROPERTIES **/

	/**
	 * Return true if this graph supports edge weights other than 1.0
	 */
	bool isWeighted() const { return weighted == Weighted::weighted; }

	/** 
	 * Return true if this graph supports directed edges.
	 */
	bool isDirected() const { return directed == Directed::directed; }

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

	void setCoordinate(node v, Point<float> value);

	Point<float>& getCoordinate(node v);

	float minCoordinate(count dim);

	float maxCoordinate(count dim);

	void initCoordinates();


	/** EDGE ATTRIBUTES **/

	/**
	 * Return edge weight.
	 *
	 * Return 0 if edge does not exist.
	 */
	edgeweight weight(node u, node v) const { return weight_impl(*this, u, v); }

	template<Directed d>
	friend edgeweight weight_impl(const BasicGraph<weighted, d>& G, node u, node v);

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
	 * Iterate over all nodes of the graph and call handler (lambda closure) as long as the condition remains true.
	 * This allows for breaking from a node loop.
	 */
	template<typename C, typename L> void forNodes(C condition, L handle) const;

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
	template<typename L> void forEdges(L handle) const { forEdges_impl(*this, handle); }

	template<Weighted w, typename L>
	void forEdges_impl(const BasicGraph<w, directed>& G, L handle);

	/**
	 * Iterate in parallel over all edges of the graph and call handler (lambda closure).
	 */
	template<typename L> void parallelForEdges(L handle) const { parallelForEdges_impl(*this, handle); }

	template<Weighted w, typename L>
	void parallelForEdges_impl(const BasicGraph<w, directed>& G, L handle);

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
	 */
	template<typename L> void forNeighborsOf(node u, L handle) const;

	template<Weighted w, typename L> 
	friend void forNeighborsOf_impl (const BasicGraph<w, directed>& G, node u, L handle);

	/**
	 * Iterate over all edge weights of a node and call handler (lamdba closure).
	 */
	template<typename L> void forWeightedNeighborsOf(node u, L handle) const;

	template<Weighted w, typename L> 
	friend void forWeightedNeighborsOf_impl (const BasicGraph<w, directed>& G, node u, L handle);

	/**
	 * Iterate over all incident edges of a node and call handler (lamdba closure).
	 */
	template<typename L> void forEdgesOf(node u, L handle) const;

	template<Weighted w, typename L> 
	friend void forEdgesOf_impl (const BasicGraph<w, directed>& G, node u, L handle);

	/**
	 * Iterate over all incident edges of a node and call handler (lamdba closure).
	 *
	 * Handle takes parameters (u, v, w) where w is the edge weight.
	 *
	 */

	template<Weighted w, typename L> 
	friend void forWeightedEdgesOf_impl (const BasicGraph<w, directed>& G, node u, L handle);

	template<typename L> void forWeightedEdgesOf(node u, L handle) const {forWeightedEdgesOf_impl(*this, u, handle);};



	/** REDUCTION ITERATORS **/
	
	/**
	 * Iterate in parallel over all nodes and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForNodes(L handle) const;

	/**
	 * Iterate in parallel over all edges and sum (reduce +) the values returned by the handler
	 */
	template<typename L> double parallelSumForWeightedEdges(L handle) const { return parallelSumForWeightedEdges(*this, handle); }

	template<Weighted w, typename L>
	friend double parallelSumForWeightedEdges_impl(const BasicGraph<w, directed>& G, L handle);


	/** GRAPH SEARCHES **/

	template<typename L> void BFSfrom(node r, L handle) const;

	template<typename L> void BFSEdgesfrom(node r, L handle) const;

	template<typename L> void DFSfrom(node r, L handle) const;

	template<typename L> void DFSEdgesfrom(node r, L handle) const;
};

} /* namespace graph_impl */

using graph_impl::BasicGraph;
// using Graph = BasicGraph<Weighted::unweighted, Directed::undirected>;
// using WeightedGraph = BasicGraph<Weighted::weighted, Directed::undirected>;
// using DirectedGraph = BasicGraph<Weighted::unweighted, Directed::directed>;
// using WeightedDirectedGraph = BasicGraph<Weighted::weighted, Directed::directed>;

template<Weighted w>
using IUndirectedGraph = BasicGraph<w, Directed::undirected>;

template<Weighted w>
using IDirectedGraph = BasicGraph<w, Directed::directed>;

template<Directed d>
using IUnweightedGraph = BasicGraph<Weighted::unweighted, d>;

template<Directed d>
using IWeigehtedGraph = BasicGraph<Weighted::weighted, d>;

} /* namespace NetworKit */

#endif /* BASICGRAPH_H_ */
