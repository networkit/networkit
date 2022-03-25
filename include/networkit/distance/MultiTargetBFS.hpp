#ifndef NETWORKIT_DISTANCE_MULTI_TARGET_BFS_HPP_
#define NETWORKIT_DISTANCE_MULTI_TARGET_BFS_HPP_

#include <networkit/distance/STSP.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Computes the shortest-path distance from a single source to multiple targets in unweighted
 * graphs.
 */
class MultiTargetBFS final : public STSP {

public:
    /**
     * Creates the MultiTargetBFS class for a graph @a G, source node @a source, and
     * multiple target nodes.
     *
     * @param G The graph.
     * @param source The source node.
     * @param targetsFirst,targetsLast Range of target nodes.
     */
    template <class InputIt>
    MultiTargetBFS(const Graph &G, node source, InputIt targetsFirst, InputIt targetsLast)
        : STSP(G, source, targetsFirst, targetsLast) {}

    void run() override;
};
} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_MULTI_TARGET_BFS_HPP_
