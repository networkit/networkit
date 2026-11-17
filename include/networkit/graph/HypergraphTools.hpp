/*
 * HypergraphTools.hpp
 *
 *  Created on: 24.05.2024
 *      Author: Fabian Brandt-Tumescheit
 */

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

} // namespace HypergraphTools

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_HYPERGRAPH_TOOLS_HPP_
