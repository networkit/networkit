/*
 * SpanningForestImpl.hpp
 *
 *  Created on: 11.09.2026
 *      Author: NetworKit contributors
 */

#ifndef NETWORKIT_GRAPH_SPANNING_FOREST_IMPL_HPP_
#define NETWORKIT_GRAPH_SPANNING_FOREST_IMPL_HPP_

#include <cassert>
#include <type_traits>
#include <vector>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/graph/BFS.hpp>

namespace NetworKit {

template <typename GraphT>
GraphT GenericSpanningForest<GraphT>::copyNodes(const GraphT &G) {
    GraphT copy(static_cast<count>(G.upperNodeIdBound()), G.isWeighted(), G.isDirected());
    for (NodeT u = 0; u < G.upperNodeIdBound(); ++u) {
        if (!G.hasNode(u)) {
            copy.removeNode(u);
        }
    }
    return copy;
}

template <typename GraphT>
void GenericSpanningForest<GraphT>::run() {
    forest = copyNodes(*G);
    std::vector<bool> visited(static_cast<index>(G->upperNodeIdBound()), false);

    const auto nodeIndex = [](NodeT u) {
        if constexpr (std::is_signed_v<NodeT>) {
            assert(u >= 0);
        }
        return static_cast<index>(u);
    };

    G->forNodes([&](NodeT source) {
        if (visited[nodeIndex(source)])
            return;

        visited[nodeIndex(source)] = true;
        Traversal::BFSEdgesFrom(*G, source, [&](NodeT u, NodeT v, EdgeWeightT weight, edgeid) {
            if (visited[nodeIndex(v)])
                return;

            visited[nodeIndex(v)] = true;
            forest.addEdge(u, v, weight);
        });
    });

    hasRun = true;
    INFO("tree edges in SpanningForest: ", forest.numberOfEdges());
}

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_SPANNING_FOREST_IMPL_HPP_
