#ifndef NETWORKIT_GRAPH_TRAVERSAL_IMPL_HPP_
#define NETWORKIT_GRAPH_TRAVERSAL_IMPL_HPP_

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace NetworKit {
namespace Traversal {
namespace Impl {

enum class Visit : uint8_t { Unvisited = 0, Visited = 1 };

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

} // namespace Impl
} // namespace Traversal
} // namespace NetworKit

#endif // NETWORKIT_GRAPH_TRAVERSAL_IMPL_HPP_
