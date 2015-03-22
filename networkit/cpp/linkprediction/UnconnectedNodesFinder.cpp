/*
 * UnconnectedNodesFinder.cpp
 *
 *  Created on: 20.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "UnconnectedNodesFinder.h"

#include <algorithm>

namespace NetworKit {

UnconnectedNodesFinder::UnconnectedNodesFinder(const Graph& G) : G(G) {
}

std::vector<std::pair<node, node>> UnconnectedNodesFinder::findAll(count k) {
  std::vector<std::pair<node, node>> missingLinks;
  G.forNodes([&](node u) {
    std::vector<std::pair<node, node>> missingAtU = findFromNode(u, k);
    for (std::pair<node, node> p : missingAtU) {
      if (std::find(missingLinks.begin(), missingLinks.end(), p) == missingLinks.end() &&
          std::find(missingLinks.begin(), missingLinks.end(), std::make_pair(p.second, p.first)) == missingLinks.end()) {
        missingLinks.push_back(p);
      }
    }
  });
  return missingLinks;
}

std::vector<std::pair<node, node>> UnconnectedNodesFinder::findFromNode(node u, count k) {
  std::vector<std::pair<node, node>> missingLinks;
  std::vector<bool> visited;
  visited.resize(G.upperNodeIdBound(), false);
  std::queue<node> q;
  q.push(u);
  visited[u] = true;
  for (count i = 1; i <= k; ++i) {
    std::queue<node> newFound;
    while (!q.empty()) {
      node u = q.front();
      q.pop();
      G.forNeighborsOf(u, [&](node v) {
        if (!visited[v]) {
          newFound.push(v);
          visited[v] = true;
        }
      });
    }
    q = newFound;
  }
  while (!q.empty()) {
    missingLinks.push_back(std::make_pair(u, q.front()));
    q.pop();
  }
  return missingLinks;
}


} // namespace NetworKit