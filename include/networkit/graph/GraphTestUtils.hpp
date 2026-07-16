#ifndef NETWORKIT_GRAPH_GRAPH_TEST_UTILS_HPP_
#define NETWORKIT_GRAPH_GRAPH_TEST_UTILS_HPP_

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit::TestUtils {

// Helper function for converting an AdjListGraph to a generic GraphT type, which is
// needed in template tests until we migrate graph readers to templates.
template <class GraphT>
GraphT ToTGraph(const Graph &G) {
    GraphT result(G.upperNodeIdBound(), G.isWeighted(), G.isDirected());
    G.forEdges([&result](node u, node v, edgeweight w) -> void {
        result.addEdge(static_cast<typename GraphT::NodeT>(u),
                       static_cast<typename GraphT::NodeT>(v),
                       static_cast<typename GraphT::EdgeWeightT>(w));
    });
    return result;
}

} // namespace NetworKit::TestUtils

#endif // NETWORKIT_GRAPH_GRAPH_TEST_UTILS_HPP_
