/*
 * BidirectionalBFS.hpp
 *
 *  Created on: 14.06.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_DISTANCE_BIDIRECTIONAL_BFS_HPP_
#define NETWORKIT_DISTANCE_BIDIRECTIONAL_BFS_HPP_

#include <queue>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/distance/STSP.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/*
 * @ingroup distance
 * Implements a bidirectional breadth-first search on a graph from two given source and target
 * nodes. Explores the graph from both the source and target nodes until the two explorations meet.
 */
class BidirectionalBFS final : public STSP {

public:
    /**
     * Creates the BidirectionalBFS class for a graph @a G, source node @a source, and
     * target node @a target.
     *
     * @param G The graph.
     * @param source The source node.
     * @param target The target node.
     * @param storePred If true, the algorithm will also store the predecessors
     * and reconstruct a shortest path from @a source and @a target.
     */
    BidirectionalBFS(const Graph &G, node source, node target, bool storePred = true)
        : STSP(G, source, target, storePred) {}

    /*
     * Perform a bidirectional BFS from the given source and target nodes.
     */
    void run() override;

    /*
     * Returns the distance (i.e., number of hops) from the source to the target node.
     *
     * @return count Number of hops from the source to the target node.
     */
    count TLX_DEPRECATED(getHops()) {
        assureFinished();
        WARN("BidirectionalBFS::getHops() is deprecated, use getDistance() instead.");
        return static_cast<count>(distance);
    }

private:
    std::vector<uint8_t> visited;
    uint8_t ts = 0;
    static constexpr uint8_t ballMask = uint8_t(1) << 7;
    std::queue<node> sQueue, sQueueNext, tQueue, tQueueNext;
};
} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_BIDIRECTIONAL_BFS_HPP_
