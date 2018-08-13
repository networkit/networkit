/*
 * BiconnectedeComponents.cpp
 *
 * Created on: March 2018
 * 		 Author: Eugenio Angriman
 */

#include <cmath>
#include <stack>

#include "BiconnectedComponents.h"

namespace NetworKit {

BiconnectedComponents::BiconnectedComponents(const Graph &G) : G(G) {
  if (G.isDirected()) {
    throw std::runtime_error(
        "Error, biconnected components cannot be computed on directed graphs.");
  }
}

void BiconnectedComponents::init() {
  n = G.numberOfNodes();
  idx = 0;
  nComp = 0;
  level.assign(n, 0);
  lowpt.assign(n, 0);
  parent.assign(n, 0);
  visited.assign(n, false);
  isRoot.assign(n, false);
  componentsOfNode.clear();
  componentsOfNode.resize(n);
}

void BiconnectedComponents::visit(node u) {
  level[u] = idx;
  lowpt[u] = idx;
  ++idx;
  visited[u] = true;
}

void BiconnectedComponents::run() {

  init();
  G.forNodes([&](node v) {
    if (visited[v]) {
      return;
    }

    isRoot[v] = true;
    std::stack<node> stack;
    std::vector<std::pair<node, node>> edgeStack;
    stack.push(v);

    while (!stack.empty()) {
      node u = stack.top();
      if (!visited[u]) {
        visit(u);
      }

      bool allVisited = true;

      for (node neighbor : G.neighbors(u)) {
        if (!visited[neighbor]) {
          allVisited = false;
          visit(neighbor);
          parent[neighbor] = u;
          stack.push(neighbor);
          edgeStack.push_back(std::make_pair(u, neighbor));
          break;
        } else if (neighbor != parent[u] && level[neighbor] < level[u]) {
          edgeStack.push_back(std::make_pair(u, neighbor));
          lowpt[u] = std::min(lowpt[u], level[neighbor]);
        }
      }

      if (allVisited) {
        stack.pop();

        if (isRoot[u]) {
          continue;
        }

        node v = parent[u];
        lowpt[v] = std::min(lowpt[v], lowpt[u]);
        if (lowpt[u] >= level[v]) {
          auto top = edgeStack.back();
          componentSizes.insert(std::make_pair(nComp, 0));
          while (level[top.first] >= level[u]) {
            edgeStack.pop_back();
            newComponent(top);
            top = edgeStack.back();
          }

          for (auto rit = edgeStack.rbegin(); rit != edgeStack.rend(); ++rit) {
            if ((*rit).first == v && (*rit).second == u) {
              newComponent(*rit);
              edgeStack.erase(std::next(rit).base());
              break;
            }
          }
          ++nComp;
        }
      }
    }
  });

  hasRun = true;
}

void BiconnectedComponents::newComponent(std::pair<node, node> e) {
  node u = e.first;
  node v = e.second;
  count prevSize = componentsOfNode[u].size() + componentsOfNode[v].size();
  componentsOfNode[u].insert(nComp);
  componentsOfNode[v].insert(nComp);
  componentSizes[nComp] +=
      componentsOfNode[u].size() + componentsOfNode[v].size() - prevSize;
}
} // namespace NetworKit
