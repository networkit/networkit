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

node randomNode(const Hypergraph &hGraph);

std::vector<node> randomNodes(const Hypergraph &hGraph, count numNodes);

edgeid randomEdge(const Hypergraph &hGraph);

std::vector<edgeid> randomEdges(const Hypergraph &hGraph, count numEdges);

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
