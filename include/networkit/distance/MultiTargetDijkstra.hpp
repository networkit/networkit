#ifndef NETWORKIT_DISTANCE_MULTI_TARGET_DIJKSTRA_HPP_
#define NETWORKIT_DISTANCE_MULTI_TARGET_DIJKSTRA_HPP_

#include <networkit/auxiliary/VectorComparator.hpp>
#include <networkit/distance/STSP.hpp>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

#include <vector>

namespace NetworKit {

/**
 * @ingroup distance
 * Computes the shortest-path distance from a single source to multiple targets in weighted graphs.
 */
class MultiTargetDijkstra final : public STSP {

public:
    /**
     * Creates the MultiTargetDijkstra class for a graph @a G, source node @a source, and
     * multiple target nodes.
     *
     * @param G The graph.
     * @param source The source node.
     * @param targetsFirst,targetsLast Range of target nodes.
     */
    template <class InputIt>
    MultiTargetDijkstra(const Graph &G, node source, InputIt targetsFirst, InputIt targetsLast)
        : STSP(G, source, targetsFirst, targetsLast) {}

    void run() override;

private:
    std::vector<edgeweight> distFromSource;
    tlx::d_ary_addressable_int_heap<node, 2, Aux::LessInVector<edgeweight>> heap{distFromSource};
};
} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_MULTI_TARGET_DIJKSTRA_HPP_
