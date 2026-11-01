#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <queue>
#include <span>
#include <stdexcept>
#include <vector>

#include <networkit/distance/DynPrunedLandmarkLabeling.hpp>
#include <networkit/distance/PrunedLandmarkLabeling.hpp>
#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

namespace {

// Returns the distances to (or from) an endpoint of a removed edge before the removal, given the
// current distances to (or from) it and to (or from) the other endpoint.
std::vector<count> distancesWithEdge(std::span<const count> dist,
                                     std::span<const count> distOtherEndpoint) {
    constexpr count infDist = std::numeric_limits<count>::max();
    std::vector<count> result(dist.begin(), dist.end());
    for (index u = 0; u < result.size(); ++u) {
        if (distOtherEndpoint[u] != infDist) {
            result[u] = std::min(result[u], distOtherEndpoint[u] + 1);
        }
    }
    return result;
}

} // namespace

DynPrunedLandmarkLabeling::DynPrunedLandmarkLabeling(const Graph &G)
    : PrunedLandmarkLabeling(G), ranks(G.upperNodeIdBound(), none) {
    for (node rank = 0; rank < nodesSortedByDegreeDesc.size(); ++rank) {
        ranks[nodesSortedByDegreeDesc[rank]] = rank;
    }
}

void DynPrunedLandmarkLabeling::update(GraphEvent e) {
    if (e.type == GraphEvent::EDGE_ADDITION) {
        addEdge(e.u, e.v);
    } else if (e.type == GraphEvent::EDGE_REMOVAL) {
        removeEdge(e.u, e.v);
    } else {
        throw std::runtime_error("Unsupported graph event " + e.toString());
    }
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
    std::ranges::fill(visited, Visit::Unvisited);
    visited[startNode] = Visit::Visited;

    std::queue<node> q0, q1;
    q0.push(startNode);

    auto visitNeighbor = [&visited = visited, &q1](node v) -> void {
        if (visited[v] == Visit::Visited) {
            return;
        }
        visited[v] = Visit::Visited;
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
    const bool directed = G->isDirected();
    const auto &labelsU = labelsOut[u];
    const auto &labelsV = directed ? labelsIn[v] : labelsOut[v];

    labelsUCopy.resize(labelsU.size());
    labelsVCopy.resize(labelsV.size());
    std::copy(labelsU.begin(), labelsU.end(), labelsUCopy.begin());
    std::copy(labelsV.begin(), labelsV.end(), labelsVCopy.begin());

    auto iterLabelsU = labelsUCopy.begin(), iterLabelsV = labelsVCopy.begin();
    const auto iterLabelsUEnd = labelsUCopy.end(), iterLabelsVEnd = labelsVCopy.end();

    // The searches must follow the order of the hubs. Otherwise, a search can reach a node with a
    // smaller rank than its root before the pair is covered, and the root becomes a hub of that
    // node. Removing edges relies on hubs always having a smaller rank than the nodes they cover.
    while (iterLabelsU != iterLabelsUEnd && iterLabelsV != iterLabelsVEnd) {
        if (iterLabelsU->node_ < iterLabelsV->node_) {
            prunedBFS(iterLabelsU->node_, v, iterLabelsU->distance_ + 1, false);
            ++iterLabelsU;
        } else if (iterLabelsU->node_ > iterLabelsV->node_) {
            prunedBFS(iterLabelsV->node_, u, iterLabelsV->distance_ + 1, directed);
            ++iterLabelsV;
        } else {
            if (directed) {
                prunedBFS(iterLabelsU->node_, v, iterLabelsU->distance_ + 1, false);
                prunedBFS(iterLabelsV->node_, u, iterLabelsV->distance_ + 1, true);
            } else if (iterLabelsU->distance_ + 1 < iterLabelsV->distance_) {
                prunedBFS(iterLabelsU->node_, v, iterLabelsU->distance_ + 1, false);
            } else {
                prunedBFS(iterLabelsV->node_, u, iterLabelsV->distance_ + 1, false);
            }

            ++iterLabelsU;
            ++iterLabelsV;
        }
    }

    while (iterLabelsU != iterLabelsUEnd) {
        prunedBFS(iterLabelsU->node_, v, iterLabelsU->distance_ + 1, false);
        ++iterLabelsU;
    }

    while (iterLabelsV != iterLabelsVEnd) {
        prunedBFS(iterLabelsV->node_, u, iterLabelsV->distance_ + 1, directed);
        ++iterLabelsV;
    }
}

std::vector<count> DynPrunedLandmarkLabeling::computeDistances(node source, bool reverse) const {
    std::vector<count> dist(G->upperNodeIdBound(), infDist);
    std::queue<node> q;
    dist[source] = 0;
    q.push(source);

    while (!q.empty()) {
        const node w = q.front();
        q.pop();

        const auto relax = [&](node z) -> void {
            if (dist[z] == infDist) {
                dist[z] = dist[w] + 1;
                q.push(z);
            }
        };

        if (reverse) {
            G->forInNeighborsOf(w, relax);
        } else {
            G->forNeighborsOf(w, relax);
        }
    }

    return dist;
}

node DynPrunedLandmarkLabeling::smallestHub(node u, node v, count dist) const {
    auto [iterLabelsU, iterLabelsUEnd] = getSourceLabelsIterators(u, G->isDirected());
    auto [iterLabelsV, iterLabelsVEnd] = getSourceLabelsIterators(v);

    while (iterLabelsU != iterLabelsUEnd && iterLabelsV != iterLabelsVEnd) {
        if (iterLabelsU->node_ < iterLabelsV->node_) {
            ++iterLabelsU;
        } else if (iterLabelsV->node_ < iterLabelsU->node_) {
            ++iterLabelsV;
        } else {
            if (iterLabelsU->distance_ + iterLabelsV->distance_ == dist) {
                return iterLabelsU->node_;
            }
            ++iterLabelsU;
            ++iterLabelsV;
        }
    }

    return none;
}

std::vector<node> DynPrunedLandmarkLabeling::detectAffectedNodes(
    node source, node target, bool reverse, std::span<const count> distToTarget,
    std::span<const count> newDistToTarget, std::span<const count> distToSource,
    std::span<Visit> isAffected) {
    // A node is visited by the search if and only if it is affected.
    std::vector<node> affected{source};
    isAffected[source] = Visit::Visited;

    const auto visitNeighbor = [&](node u) -> void {
        if (isAffected[u] == Visit::Visited) {
            return;
        }

        bool isAffectedU = newDistToTarget[u] != distToTarget[u];
        if (!isAffectedU) {
            const node rank = reverse ? smallestHub(u, target, distToTarget[u])
                                      : smallestHub(target, u, distToTarget[u]);
            assert(rank != none);
            const node hub = nodesSortedByDegreeDesc[rank];
            isAffectedU =
                isAffected[hub] == Visit::Visited
                || ((hub == u || hub == target) && distToTarget[u] == distToSource[u] + 1);
        }

        if (isAffectedU) {
            isAffected[u] = Visit::Visited;
            affected.push_back(u);
        }
    };

    // The affected vector doubles as the queue of the BFS.
    for (index i = 0; i < affected.size(); ++i) {
        if (reverse) {
            G->forInNeighborsOf(affected[i], visitNeighbor);
        } else {
            G->forNeighborsOf(affected[i], visitNeighbor);
        }
    }

    return affected;
}

void DynPrunedLandmarkLabeling::removeHubs(std::vector<std::vector<Label>> &labels,
                                           std::span<const node> nodes,
                                           std::span<const Visit> isHub) {
    for (const node u : nodes) {
        std::erase_if(labels[u], [&](const Label &label) -> bool {
            const node hub = nodesSortedByDegreeDesc[label.node_];
            return hub != u && isHub[hub] == Visit::Visited;
        });
    }
}

void DynPrunedLandmarkLabeling::insertLabel(std::vector<Label> &labelsU, node hub, count dist) {
    labelsU.emplace(std::ranges::lower_bound(labelsU, hub, {}, &Label::node_), hub, dist);
}

void DynPrunedLandmarkLabeling::restoreCover(node root, std::span<const Visit> isTarget,
                                             bool reverse, bool ordered) {
    // For directed graphs, labelsIn stores the distances to the hubs and labelsOut the distances
    // from the hubs.
    auto &labelsToHub = G->isDirected() ? labelsIn : labelsOut;

    std::ranges::fill(visited, Visit::Unvisited);
    visited[root] = Visit::Visited;

    std::queue<node> q0, q1;
    q0.push(root);
    count level = 0;

    const auto visitNeighbor = [&](node w) -> void {
        if (visited[w] == Visit::Unvisited) {
            visited[w] = Visit::Visited;
            q1.push(w);
        }
    };

    do {
        do {
            const node v = q0.front();
            q0.pop();

            if (ordered && ranks[v] < ranks[root]) {
                continue;
            }

            if (isTarget[v] == Visit::Visited
                && level < (reverse ? queryImpl(v, root) : queryImpl(root, v))) {
                if (ranks[v] < ranks[root]) {
                    insertLabel(reverse ? labelsOut[root] : labelsToHub[root], ranks[v], level);
                } else {
                    insertLabel(reverse ? labelsToHub[v] : labelsOut[v], ranks[root], level);
                }
            }

            if (reverse) {
                G->forInNeighborsOf(v, visitNeighbor);
            } else {
                G->forNeighborsOf(v, visitNeighbor);
            }
        } while (!q0.empty());

        ++level;
        std::swap(q0, q1);
    } while (!q0.empty());
}

void DynPrunedLandmarkLabeling::removeEdge(node x, node y) {
    // Self-loops and parallel edges do not affect distances.
    if (x == y || G->hasEdge(x, y)) {
        return;
    }

    const bool directed = G->isDirected();

    // The distances before the removal follow from the current ones: the shortest paths that used
    // the edge started with it (if from x) or ended with it (if to y). In undirected graphs, the
    // distances from a node are the distances to it.
    const auto newToX = computeDistances(x, directed);
    const auto newToY = computeDistances(y, directed);
    const auto newFromXDirected = directed ? computeDistances(x, false) : std::vector<count>{};
    const auto newFromYDirected = directed ? computeDistances(y, false) : std::vector<count>{};
    const std::span<const count> newFromX = directed ? newFromXDirected : newToX;
    const std::span<const count> newFromY = directed ? newFromYDirected : newToY;

    const auto toY = distancesWithEdge(newToY, newToX);
    const auto fromX = distancesWithEdge(newFromX, newFromY);
    // In directed graphs, no shortest path to x or from y can use the edge (x, y).
    const std::span<const count> toX = directed ? newToX : fromX;
    const std::span<const count> fromY = directed ? newFromYDirected : toY;

    // Phase 1: nodes with a shortest path through the removed edge that is induced by the labels.
    std::vector<Visit> isAffectedX(G->upperNodeIdBound(), Visit::Unvisited);
    std::vector<Visit> isAffectedY(G->upperNodeIdBound(), Visit::Unvisited);
    const auto affectedX = detectAffectedNodes(x, y, directed, toY, newToY, toX, isAffectedX);
    const auto affectedY = detectAffectedNodes(y, x, false, fromX, newFromX, fromY, isAffectedY);

    // Phase 2: only labels between nodes of the two sets can be outdated.
    removeHubs(directed ? labelsIn : labelsOut, affectedX, isAffectedY);
    removeHubs(labelsOut, affectedY, isAffectedX);

    // If the endpoints are disconnected, so are all pairs of nodes of the two sets. This does not
    // hold for directed graphs, where the two sets can be connected by paths that avoid x and y.
    if (!directed && newToY[x] == infDist) {
        return;
    }

    // Phase 3: Order-Restore is faster unless both sets are large (roughly n / log n according to
    // the experiments in the paper).
    const count n = G->numberOfNodes();
    const count minAffected = std::min(affectedX.size(), affectedY.size());
    if (static_cast<double>(minAffected) <= static_cast<double>(n) / std::log2(n)) {
        std::vector<node> roots(affectedX);
        roots.insert(roots.end(), affectedY.begin(), affectedY.end());
        std::ranges::sort(roots, {}, [this](node u) -> node { return ranks[u]; });
        roots.erase(std::unique(roots.begin(), roots.end()), roots.end());

        for (const node root : roots) {
            if (isAffectedX[root] == Visit::Visited) {
                restoreCover(root, isAffectedY, /*reverse=*/false, /*ordered=*/true);
            }
            if (isAffectedY[root] == Visit::Visited) {
                restoreCover(root, isAffectedX, /*reverse=*/directed, /*ordered=*/true);
            }
        }
    } else if (affectedX.size() <= affectedY.size()) {
        for (const node root : affectedX) {
            restoreCover(root, isAffectedY, /*reverse=*/false, /*ordered=*/false);
        }
    } else {
        for (const node root : affectedY) {
            restoreCover(root, isAffectedX, /*reverse=*/directed, /*ordered=*/false);
        }
    }
}

} // namespace NetworKit
