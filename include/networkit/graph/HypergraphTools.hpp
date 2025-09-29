#ifndef NETWORKIT_GRAPH_HYPERGRAPH_TOOLS_HPP_
#define NETWORKIT_GRAPH_HYPERGRAPH_TOOLS_HPP_

#include <networkit/Globals.hpp>
#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {
namespace HypergraphTools {

// /**
//  * @brief Computes the linegraph of a hypergraph.
//  *
//  * @param hGraph Input hypergraph.
//  * @return Graph Linegraph.
//  */
// Graph computeLineGraph(const Hypergraph &hGraph);

/**
 * Returns a uniformly at random chosen node from the hypergraph.
 *
 * @param hGraph The Hypergraph.
 * @return node Randomly chosen node.
 */
node randomNode(const Hypergraph &hGraph);

/**
 * Returns unique uniformly at random chosen nodes from the Hypergraph.
 *
 * @param hGraph The Hypergraph.
 * @param numNodes Number of nodes, which should be returned.
 * @return std::vector<node> Randomly chosen nodes.
 */
std::vector<node> randomNodes(const Hypergraph &hGraph, count numNodes);

/**
 * Returns a uniformly at random chosen edge from the hypergraph.
 *
 * @param hGraph The Hypergraph.
 * @return edge Randomly chosen edge.
 */
edgeid randomEdge(const Hypergraph &hGraph);

/**
 * Returns unique uniformly at random chosen edges from the Hypergraph.
 *
 * @param hGraph The Hypergraph.
 * @param numEdges Number of edges, which should be returned.
 * @return std::vector<edgeid> Randomly chosen edges.
 */
std::vector<edgeid> randomEdges(const Hypergraph &hGraph, count numEdges);

/**
 * Computes the max edge order (cardinality) for a given Hypergraph.
 *
 * @param hGraph The Hypergraph.
 * @return Maximum edge order.
 */
count maxEdgeOrder(const Hypergraph &hGraph);

/**
 * Computes the maximum unweighted degree for a given Hypergraph.
 *
 * @param hGraph The Hypergraph.
 * @return Maximum unweighted node degree.
 */
count maxDegree(const Hypergraph &hGraph);

/**
 * Computes the maximum weighted degree for a given Hypergraph.
 *
 * @param hGraph The Hypergraph.
 * @return Maximum weighted node degree.
 */
edgeweight maxWeightedDegree(const Hypergraph &hGraph);

/**
 * Returns the intersection of two hyperedges as a set
 *
 * @param hypergraph The Hypergraph.
 * @param eid1 first hyperedge
 * @param eid2 second hyperedge
 * @return intersection
 */
std::unordered_set<node> getIntersection(Hypergraph &hGraph, edgeid eid1, edgeid eid2);

/**
 * Returns the size of the intersection of two hyperedges
 *
 * @param hypergraph The Hypergraph.
 * @param eid1 first hyperedge
 * @param eid2 second hyperedge
 * @return intersection size
 */
count getIntersectionSize(Hypergraph &hGraph, edgeid eid1, edgeid eid2);

/**
 * Converts a hypergraph into its clique expansion (a simple graph)
 *
 * @param hGraph The Hypergraph.
 * @return cliqueExpansion The clique expansion.
 */
Graph cliqueExpansion(Hypergraph &hGraph);

/**
 * Converts a hypergraph into its line expansion (an attributed graph).
 * The attributes in the graph contain the original node and edge id from the Hypergraph.
 *
 * @param hGraph The Hypergraph.
 * @return lineExpansion The line expansion.
 */
Graph lineExpansion(Hypergraph &hGraph);

/**
 * Converts a hypergraph into its line graph (a simple graph) with optional weights.
 * The weights are computed based on: 1/3 * (union + union/intersection) of two hyperedges.
 * The formula is taken from "Distances in Higher-Order Networks and the Metric Structure of
 * Hypergraphs" (E. Vasilyeva et al.) Link: https://doi.org/10.3390/e25060923
 *
 * @param hGraph The Hypergraph.
 * @param weighted weighted Optional: If set to true, the edges of the line graph are weighted based
 * on their order and intersection size. Default value: false.
 * @return lineGraph The line graph.
 */
Graph lineGraph(Hypergraph &hGraph, bool weighted = false);

// template <typename Matrix>
// Matrix computeIncidenceMatrix(const Hypergraph &hGraph);

// template <>
// CSRMatrix computeIncidenceMatrix(const Hypergraph &hGraph);

// template <>
// DenseMatrix computeIncidenceMatrix(const Hypergraph &hGraph);

// template <>
// DynamicMatrix computeIncidenceMatrix(const Hypergraph &hGraph);

} // namespace HypergraphTools

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_HYPERGRAPH_TOOLS_HPP_
