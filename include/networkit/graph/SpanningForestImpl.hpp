/*
 * SpanningForestImpl.hpp
 *
 *  Created on: 11.09.2026
 *      Author: NetworKit contributors
 */

#ifndef NETWORKIT_GRAPH_SPANNING_FOREST_IMPL_HPP_
#define NETWORKIT_GRAPH_SPANNING_FOREST_IMPL_HPP_

#include <cassert>
#include <queue>
#include <type_traits>
#include <vector>

#include <networkit/auxiliary/Log.hpp>

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

        std::queue<NodeT> queue;
        queue.push(source);
        visited[nodeIndex(source)] = true;

        while (!queue.empty()) {
            const NodeT u = queue.front();
            queue.pop();

            if (G->isWeighted()) {
                for (const auto &[v, weight] : G->weightNeighborRange(u)) {
                    if (visited[nodeIndex(v)])
                        continue;

                    visited[nodeIndex(v)] = true;
                    forest.addEdge(u, v, weight);
                    queue.push(v);
                }
            } else {
                for (const NodeT v : G->neighborRange(u)) {
                    if (visited[nodeIndex(v)])
                        continue;

                    visited[nodeIndex(v)] = true;
                    forest.addEdge(u, v);
                    queue.push(v);
                }
            }
        }
    });

    hasRun = true;
    INFO("tree edges in SpanningForest: ", forest.numberOfEdges());
}

} // namespace NetworKit

#endif // NETWORKIT_GRAPH_SPANNING_FOREST_IMPL_HPP_
