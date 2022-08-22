#include <networkit/distance/BFS.hpp>
#include <networkit/distance/MultiTargetBFS.hpp>

#include <limits>
#include <queue>
#include <unordered_set>

namespace NetworKit {

void MultiTargetBFS::run() {
    distances.clear();
    distances.reserve(targets.size());
    targetIdx.clear();

    std::unordered_set<node> targetsSet(targets.begin(), targets.end());
    count level = 0;

    BFS bfs(*G, source);
    bfs.run();
    std::queue<node> frontier, next;
    frontier.push(source);

    std::vector<bool> visited(G->upperNodeIdBound());

    auto visitNode = [&](node u) -> bool {
        if (visited[u]) // Already visited, nothing to do
            return false;

        // Check if we are visiting a target node
        const auto it = targetsSet.find(u);
        if (it != targetsSet.end()) {
            targetIdx.emplace(u, distances.size());
            distances.push_back(static_cast<edgeweight>(level));
            targetsSet.erase(it);
            if (targetsSet.empty())
                return true; // All target nodes visited, stop exploration
        }

        visited[u] = true;
        next.push(u);
        return false;
    };

    auto toNextLevel = [&]() -> void {
        ++level;
        std::swap(frontier, next);
    };

    visitNode(source);
    toNextLevel();

    while (!targetsSet.empty() && !frontier.empty()) {
        do {
            const auto front = frontier.front();
            frontier.pop();

            for (node neighbor : G->neighborRange(front))
                if (visitNode(neighbor))
                    break;
        } while (!frontier.empty());

        toNextLevel();
    }

    // Mark targets not reached as unreachable
    for (const node target : targetsSet) {
        targetIdx.emplace(target, distances.size());
        distances.emplace_back(std::numeric_limits<edgeweight>::max());
    }

    hasRun = true;
}
} // namespace NetworKit
