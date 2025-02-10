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
    isGraphBiPartite = true;
    hasRun = true;
    graph->forNodesWhile([&]() { return isGraphBiPartite == true; },
                         [&](node u) {
                             if (colors[u] != uncolored)
                                 return;
                             colors[u] = black;
                             std::queue<node> queue;
                             queue.push(u);
                             do {
                                 const node currentNode = queue.front();
                                 queue.pop();
                                 for (const auto neighbor : graph->neighborRange(currentNode)) {
                                     if (colors[neighbor] == uncolored) {
                                         colors[neighbor] = red - colors[currentNode];
                                         queue.push(neighbor);
                                     } else if (colors[neighbor] == colors[currentNode]) {
                                         isGraphBiPartite = false;
                                         return;
                                     }
                                 }
                             } while (!queue.empty());
                         });
}

} // namespace NetworKit
