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
    template <class InputIt>
    MultiTargetBFS(const Graph &G, node source, InputIt targetsFirst, InputIt targetsLast)
        : STSP(G, source, targetsFirst, targetsLast) {}

    void run() override;
};
} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_MULTI_TARGET_BFS_HPP_
