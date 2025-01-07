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
#include <networkit/matching/BSuitorMatcher.hpp>

namespace NetworKit {
class DynamicBSuitorMatcher final : public BSuitorMatcher {
    enum class Operation { Insert, Remove };

    bool isBetterMatch(node u, node v, edgeweight ew) const noexcept {
        const auto currentMatch = suitors[u].min;
        bool isBetterMatch = currentMatch.id == none || currentMatch.weight < ew
                             || (currentMatch.weight == ew && v < currentMatch.id);
        return isBetterMatch;
    }

    void trackUpdatePath(size_t batchId, node start, bool recursiveCall = false);
    void processEdgeInsertion(const WeightedEdge &edge);
    void processEdgeRemoval(const Edge &edge);

public:
    DynamicBSuitorMatcher(const Graph &G, const std::vector<count> &b) : BSuitorMatcher(G, b) {}
    DynamicBSuitorMatcher(const Graph &G, count b = 1) : BSuitorMatcher(G, b) {}
    DynamicBSuitorMatcher(const Graph &G, const std::string &path) : BSuitorMatcher(G, 1) {}

    void addEdge(WeightedEdge &edge) {
        std::vector<WeightedEdge> edges{edge};
        addEdges(edges);
    }
    void addEdges(std::vector<WeightedEdge> &edges, bool sort = true);

    void removeEdge(Edge &edge) {
        std::vector<Edge> edges{edge};
        removeEdges(edges);
    }
    void removeEdges(std::vector<Edge> &edges);
};

} // namespace NetworKit
#endif // NETWORKIT_MATCHING_DYNAMIC_B_SUITOR_MATCHER_HPP_
