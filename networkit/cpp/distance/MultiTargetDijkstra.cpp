#include <networkit/distance/MultiTargetDijkstra.hpp>

#include <limits>
#include <unordered_set>

namespace NetworKit {

void MultiTargetDijkstra::run() {
    distances.clear();
    distances.reserve(targets.size());
    targetIdx.clear();
    heap.clear();

    const auto infDist = std::numeric_limits<edgeweight>::max();
    distFromSource.assign(G->upperNodeIdBound(), infDist);
    distFromSource[source] = 0;

    heap.reserve(G->upperNodeIdBound());
    heap.push(source);

    std::unordered_set<node> targetsSet(targets.begin(), targets.end());

    do {
        const auto top = heap.extract_top();

        const auto it = targetsSet.find(top);
        if (it != targetsSet.end()) {
            targetIdx.emplace(top, distances.size());
            distances.push_back(distFromSource[top]);
            targetsSet.erase(it);
            if (targetsSet.empty())
                break;
        }

        G->forNeighborsOf(top, [&](node u, edgeweight weight) {
            const auto distU = distFromSource[top] + weight;
            if (distFromSource[u] > distU) {
                distFromSource[u] = distU;
                heap.update(u);
            }
        });
    } while (!heap.empty());

    // Mark targets not reached as unreachable
    for (const node target : targetsSet) {
        targetIdx.emplace(target, distances.size());
        distances.emplace_back(std::numeric_limits<edgeweight>::max());
    }

    hasRun = true;
}

} // namespace NetworKit
