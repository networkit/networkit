/*
 * MissingLinksFinder.cpp
 *
 *  Created on: 20.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "MissingLinksFinder.h"

#include <algorithm>
#include <random>

namespace NetworKit {

MissingLinksFinder::MissingLinksFinder(const Graph& G) : G(G) {
}

std::vector<std::pair<node, node>> MissingLinksFinder::findRandomly(count k, count limit) {
  std::set<std::pair<node, node>> missingLinks;
  std::vector<node> nodes = G.nodes();
  std::random_device randomDevice;
  std::mt19937 generator(randomDevice());
  count numMissingLinks = 0;
  // Determine the number of missing links with distance k
  /*#pragma omp parallel for
  for (index i = 0; i < nodes.size(); ++i) {
    std::vector<std::pair<node, node>> missingAtU = findFromNode(nodes[i], k);
    #pragma omp atomic
    numMissingLinks += missingAtU.size();
  }
  // Every edge will be counted twice, so adjust the number accordingly
  numMissingLinks /= 2;
  if (limit > numMissingLinks) {
    return findAll(k);
  }*/
  while (missingLinks.size() < limit) {
    std::vector<std::pair<node, node>> missingAtU;
    // Make sure to actually get a node that has missing links with distance k
    while (missingAtU.size() == 0) {
      node randomNode = G.randomNode();
      missingAtU = findFromNode(randomNode, k);
      INFO("Random node = ", randomNode);
    }
    std::uniform_int_distribution<> distribution(0, missingAtU.size() - 1);
    std::pair<node, node> edge = missingAtU[distribution(generator)];
    missingLinks.insert(edge.first < edge.second ? edge : std::make_pair(edge.second, edge.first));
  }
  std::vector<std::pair<node, node>> result(limit);
  result.assign(missingLinks.begin(), missingLinks.end());
  return result;
}

std::vector<std::pair<node, node>> MissingLinksFinder::findAll(count k) {
  std::vector<std::pair<node, node>> missingLinks;
  std::vector<node> nodes = G.nodes();
  #pragma omp parallel
  {
    std::vector<std::pair<node, node>> missingLinksPrivate;
    #pragma omp for nowait
    for (index i = 0; i < nodes.size(); ++i) {
      std::vector<std::pair<node, node>> missingAtU = findFromNode(nodes[i], k);
      // Discard all node-pairs of the form u > v. This removes all duplicates that result from undirected edges.
      missingAtU.erase(std::remove_if(std::begin(missingAtU), std::end(missingAtU),
          [&](std::pair<node, node> p) { return p.first >= p.second; }), std::end(missingAtU));
      missingLinksPrivate.insert(missingLinksPrivate.end(), missingAtU.begin(), missingAtU.end());
    }
    #pragma omp critical
    missingLinks.insert(missingLinks.end(), missingLinksPrivate.begin(), missingLinksPrivate.end());
  }
  DEBUG("Found ", missingLinks.size(), " missing links with distance ", k, ".");
  std::sort(missingLinks.begin(), missingLinks.end());
  return missingLinks;
}

std::vector<std::pair<node, node>> MissingLinksFinder::findFromNode(node u, count k) {
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