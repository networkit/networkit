/*
 * DynamicBSuitorMatcher.hpp
 *
 *  Created on: 06.01.2025
 *      Author: Fabian Brandt-Tumescheit
 *              Frieda Gerharz
 */

#ifndef NETWORKIT_MATCHING_DYNAMIC_B_SUITOR_MATCHER_HPP_
#define NETWORKIT_MATCHING_DYNAMIC_B_SUITOR_MATCHER_HPP_

#include <set>
#include <unordered_map>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/matching/BSuitorMatcher.hpp>

namespace NetworKit {
class DynamicBSuitorMatcher final : public BSuitorMatcher, public DynAlgorithm {
    enum class Operation { Insert, Remove };

public:
    /**
     * Implementation from the algorithm from "A Fully-dynamic Approximation Algorithm for Maximum
     * Weight b-Matchings in Graphs" from Proceedings of The Thirteenth International Conference on
     * Complex Networks and their Applications 2024 by Fabian Brandt-Tumescheit, Frieda Gerharz and
     * Henning Meyerhenke. The algorithm dynamically updates the b-matching based on the b-Suitor
     * algorithm by Khan et al. The solution is the same as a complete recomputation of the b-Suitor
     * b-matching.
     *
     * @param   G	The graph,
     * @param   b   List of b-values for each node in the graph.
     */
    DynamicBSuitorMatcher(const Graph &G, const std::vector<count> &b) : BSuitorMatcher(G, b) {}

    /**
     * @param   G	The graph,
     * @param   b   Set the same b-value for all nodes. Optional and defaults to 1.
     */
    DynamicBSuitorMatcher(const Graph &G, count b = 1) : BSuitorMatcher(G, b) {}

    /**
     * Updates the b-matching after an edge insertion or deletion on the graph.
     * Notice: Supported events include edge insertion and deletion.
     *
     * @param e The update event.
     */
    void update(GraphEvent e) override;

    /**
     * Updates the b-matching after a batch of edge insertions or deletions on the graph.
     * Notice: Supported events include edge insertion and deletion.
     *
     * @param e The batch of update events.
     */
    void updateBatch(const std::vector<GraphEvent> &batch) override;

private:
    // helper function
    bool isBetterMatch(node u, node v, edgeweight ew) const noexcept {
        const MatchingNode currentMatch = suitors[u].min;
        return currentMatch.id == none || currentMatch.weight < ew
               || (currentMatch.weight == ew && v < currentMatch.id);
    }

    void addEdge(const GraphEvent &event);

    void processEdgeInsertion(const GraphEvent &event);

    void removeEdge(const GraphEvent &event);

    void processEdgeRemoval(const GraphEvent &event);

    void trackUpdatePath(node start);
};

} // namespace NetworKit
#endif // NETWORKIT_MATCHING_DYNAMIC_B_SUITOR_MATCHER_HPP_
