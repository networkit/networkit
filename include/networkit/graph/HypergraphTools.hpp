#ifndef NETWORKIT_GRAPH_HYPERGRAPH_TOOLS_HPP_
#define NETWORKIT_GRAPH_HYPERGRAPH_TOOLS_HPP_

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {
namespace HypergraphTools {

/**
 * @brief Computes the linegraph of a hypergraph.
 *
 * @param hGraph Input hypergraph.
 * @return Graph Linegraph.
 */
Graph computeLineGraph(const Hypergraph &hGraph);

template <typename Matrix>
Matrix computeIncidenceMatrix(const Hypergraph &hGraph);

template <>
CSRMatrix computeIncidenceMatrix(const Hypergraph &hGraph);

template <>
DenseMatrix computeIncidenceMatrix(const Hypergraph &hGraph);

template <>
DynamicMatrix computeIncidenceMatrix(const Hypergraph &hGraph);

} // namespace HypergraphTools

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_GRAPH_TOOLS_HPP_
