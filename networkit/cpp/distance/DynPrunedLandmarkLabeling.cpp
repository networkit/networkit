#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <utility>

#include <networkit/distance/DynPrunedLandmarkLabeling.hpp>
#include <networkit/distance/PrunedLandmarkLabeling.hpp>
#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

void DynPrunedLandmarkLabeling::update(GraphEvent e) {
    if (e.type == GraphEvent::EDGE_ADDITION)
        addEdge(e.u, e.v);
    else
        throw std::runtime_error("Unsupported graph event " + e.toString());
}

void DynPrunedLandmarkLabeling::sortUpdatedLabels(bool reverse) {
    for (const node u : updatedNodes) {
        auto &labelsU = reverse ? labelsIn[u] : labelsOut[u];
        if (labelsU.size() < 2) // Nothing to be sorted
            continue;

        // Insertion sort step: insert the new label in the right position within the vector.
        auto lb =
            std::lower_bound(labelsU.begin(), labelsU.end() - 2, labelsU.back(),
                             [](const auto &l1, const auto &l2) { return l1.node_ < l2.node_; });

        if (lb->node_ == labelsU.back().node_ && lb->distance_ > labelsU.back().distance_) {
            // Overwrite label
            *lb = labelsU.back();
            labelsU.pop_back();
        } else if (lb->node_ > labelsU.back().node_) {
            // Insert new label
            index i = labelsU.size() - 1, start = (lb - labelsU.begin());
            const auto newLabel = labelsU.back();
            while (i > start) {
                labelsU[i] = labelsU[i - 1];
                --i;
            }
            labelsU[i] = newLabel;
        }
    }
}

void DynPrunedLandmarkLabeling::prunedBFS(node k, node startNode, count bfsLevel, bool reverse) {
    const node root = nodesSortedByDegreeDesc[k];

    updatedNodes.clear();
    std::fill(visited.begin(), visited.end(), false);
    visited[startNode] = true;

    std::queue<node> q0, q1;
    q0.push(startNode);

    auto visitNeighbor = [&visited = visited, &q1](node v) -> void {
        if (visited[v])
            return;
        visited[v] = true;
        q1.push(v);
    };

    do {
        do {
            const node u = q0.front();
            q0.pop();

            if (reverse) {
                if (queryImpl(u, root, k) <= bfsLevel)
                    continue;
            } else {
                if (queryImpl(root, u, k) <= bfsLevel)
                    continue;
            }

            updatedNodes.push_back(u);

            if (reverse) {
                labelsIn[u].emplace_back(k, bfsLevel);
                G->forInNeighborsOf(u, visitNeighbor);
            } else {
                labelsOut[u].emplace_back(k, bfsLevel);
                G->forNeighborsOf(u, visitNeighbor);
            }
        } while (!q0.empty());

        ++bfsLevel;
        std::swap(q0, q1);
    } while (!q0.empty());

    sortUpdatedLabels(reverse);
}

void DynPrunedLandmarkLabeling::addEdge(node u, node v) {
    const auto &labelsU = labelsOut[u];
    const auto &labelsV = G->isDirected() ? labelsIn[v] : labelsOut[v];

    labelsUCopy.resize(labelsU.size());
    labelsVCopy.resize(labelsV.size());
    std::copy(labelsU.begin(), labelsU.end(), labelsUCopy.begin());
    std::copy(labelsV.begin(), labelsV.end(), labelsVCopy.begin());

    auto iterLabelsU = labelsUCopy.begin(), iterLabelsV = labelsVCopy.begin();
    const auto iterLabelsUEnd = labelsUCopy.end(), iterLabelsVEnd = labelsVCopy.end();

    if (!G->isDirected()) {
        do {
            if (iterLabelsU->node_ < iterLabelsV->node_) {
                prunedBFS(iterLabelsU->node_, v, iterLabelsU->distance_ + 1, /*reverse=*/false);
                ++iterLabelsU;
            } else if (iterLabelsU->node_ > iterLabelsV->node_) {
                prunedBFS(iterLabelsV->node_, u, iterLabelsV->distance_ + 1, /*reverse=*/false);
                ++iterLabelsV;
            } else {
                if (iterLabelsU->distance_ + 1 < iterLabelsV->distance_)
                    prunedBFS(iterLabelsU->node_, v, iterLabelsU->distance_ + 1, /*reverse=*/false);
                else
                    prunedBFS(iterLabelsV->node_, u, iterLabelsV->distance_ + 1, /*reverse=*/false);

                ++iterLabelsU;
                ++iterLabelsV;
            }
        } while (iterLabelsU != iterLabelsUEnd && iterLabelsV != iterLabelsVEnd);
    }

    while (iterLabelsU != iterLabelsUEnd) {
        prunedBFS(iterLabelsU->node_, v, iterLabelsU->distance_ + 1, /*reverse=*/false);
        ++iterLabelsU;
    }

    while (iterLabelsV != iterLabelsVEnd) {
        prunedBFS(iterLabelsV->node_, u, iterLabelsV->distance_ + 1, /*reverse=*/true);
        ++iterLabelsV;
    }
}

} // namespace NetworKit
