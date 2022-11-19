#include <algorithm>
#include <queue>
#include <vector>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/distance/PrunedLandmarkLabeling.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

PrunedLandmarkLabeling::PrunedLandmarkLabeling(const Graph &G)
    : G(&G), nodesSortedByDegreeDesc(G.nodeRange().begin(), G.nodeRange().end()) {

    if (G.isWeighted())
        WARN("This algorithm ignores edge weights.");

    if (G.isDirected())
        Aux::Parallel::sort(nodesSortedByDegreeDesc.begin(), nodesSortedByDegreeDesc.end(),
                            [&G](node u, node v) {
                                count degU = G.degree(u), degV = G.degree(v);
                                if (degU != degV)
                                    return degU > degV;
                                return G.degreeIn(u) > G.degreeIn(v);
                            });
    else
        Aux::Parallel::sort(nodesSortedByDegreeDesc.begin(), nodesSortedByDegreeDesc.end(),
                            [&G](node u, node v) { return G.degree(u) > G.degree(v); });

    visited.resize(G.upperNodeIdBound());

    labelsOut.resize(G.upperNodeIdBound());
    if (G.isDirected())
        labelsIn.resize(G.upperNodeIdBound());

    labelsUCopy.reserve(G.upperNodeIdBound());
    labelsVCopy.reserve(G.upperNodeIdBound());
}

template <bool Reverse>
void PrunedLandmarkLabeling::prunedBFS(node root, node rankOfRootNode) {
    std::fill(visited.begin(), visited.end(), false);
    visited[root] = true;

    std::queue<node> q1, q2;
    q1.push(root);

    index level = 0;

    const auto visitNeighbor = [&](node v) -> void {
        if (visited[v])
            return;
        visited[v] = true;
        q2.push(v);
    };

    do {
        do {
            const node u = q1.front();
            q1.pop();
            if constexpr (Reverse) {
                if (u != root && queryImpl(u, root) <= level)
                    continue;
            } else {
                if (u != root && queryImpl(root, u) <= level)
                    continue;
            }

            if constexpr (Reverse) {
                labelsIn[u].emplace_back(rankOfRootNode, level);
                G->forInNeighborsOf(u, visitNeighbor);
            } else {
                labelsOut[u].emplace_back(rankOfRootNode, level);
                G->forNeighborsOf(u, visitNeighbor);
            }

        } while (!q1.empty());

        ++level;
        std::swap(q1, q2);
    } while (!q1.empty());
}

void PrunedLandmarkLabeling::run() {
    index rankOfRootNode = 0;
    for (node root : nodesSortedByDegreeDesc) {
        prunedBFS(root, rankOfRootNode);
        if (G->isDirected())
            prunedBFS<true>(root, rankOfRootNode);
        ++rankOfRootNode;
    }

    hasRun = true;
}

count PrunedLandmarkLabeling::queryImpl(node u, node v) const {
    if (u == v)
        return 0;

    auto [iterLabelsU, iterLabelsUEnd] = getSourceLabelsIterators(u, G->isDirected());
    auto [iterLabelsV, iterLabelsVEnd] = getSourceLabelsIterators(v);

    count result = infDist;

    while (iterLabelsU != iterLabelsUEnd && iterLabelsV != iterLabelsVEnd) {
        if (iterLabelsU->node_ < iterLabelsV->node_)
            ++iterLabelsU;
        else if (iterLabelsV->node_ < iterLabelsU->node_)
            ++iterLabelsV;
        else {
            result = std::min(result, iterLabelsU->distance_ + iterLabelsV->distance_);

            ++iterLabelsU;
            ++iterLabelsV;
        }
    }

    return result;
}

count PrunedLandmarkLabeling::query(node u, node v) const {
    assureFinished();
    return queryImpl(u, v);
}

} // namespace NetworKit
