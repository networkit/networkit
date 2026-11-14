/*  TraversalImpl.hpp
 *
 *  Created on: 25.10.2026
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */

#ifndef NETWORKIT_GRAPH_TRAVERSAL_IMPL_HPP_
#define NETWORKIT_GRAPH_TRAVERSAL_IMPL_HPP_

#include <networkit/graph/GraphConcepts.hpp>

#include <cassert>
#include <cstddef>
#include <type_traits>

namespace NetworKit::Traversal::Impl {

template <class F, class NodeT, class DistanceT>
concept TraversalHandle =
    GraphNode<NodeT> && (std::invocable<F &, NodeT> || std::invocable<F &, NodeT, DistanceT>);

template <typename NodeT>
std::size_t nodeIndex(NodeT u) {
    if constexpr (std::is_signed_v<NodeT>) {
        assert(u >= 0);
    }
    return static_cast<std::size_t>(u);
}

template <class GraphT, typename NodeT>
bool isValidNode(const GraphT &G, NodeT u) {
    if constexpr (std::is_signed_v<NodeT>) {
        if (u < 0)
            return false;
    }
    return G.hasNode(u);
}

} // namespace NetworKit::Traversal::Impl

#endif // NETWORKIT_GRAPH_TRAVERSAL_IMPL_HPP_
