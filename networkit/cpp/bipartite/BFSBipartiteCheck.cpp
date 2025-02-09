/*  BFSBipartiteCheck.cpp
 *
 *	Created on: 09.02.2025
 *  Authors: Andreas Scharf (andreas.b.scharf@gmail.com)
 *
 */
#include <queue>
#include <networkit/bipartite/BFSBipartiteCheck.hpp>
namespace NetworKit {

void BFSBipartiteCheck::run() {
    constexpr signed char uncolored{-1};
    constexpr signed char black{};
    constexpr signed char red{1};
    std::vector<signed char> colors(graph->numberOfNodes(), uncolored);
    std::queue<node> q;
    isGraphBiPartite = true;
    for (node u = 0; u < graph->numberOfNodes(); u++) {
        if (colors[u] == uncolored) {
            colors[u] = black;
            q.push(u);
            do {
                const node currentNode = q.front();
                q.pop();
                for (const auto neighbor : graph->neighborRange(currentNode)) {
                    if (colors[neighbor] == uncolored) {
                        colors[neighbor] = red - colors[currentNode];
                        q.push(neighbor);
                    } else if (colors[neighbor] == colors[currentNode]) {
                        isGraphBiPartite = false;
                        return;
                    }
                }
            } while (!q.empty());
        }
    }
}

} // namespace NetworKit
