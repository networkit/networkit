/*
 * IncompleteDijkstra.cpp
 *
 *  Created on: 15.07.2014
 *      Author: dhoske
 */

#include <limits>
#include <type_traits>

#include "IncompleteDijkstra.h"

using namespace std;

namespace NetworKit {

void IncompleteDijkstra::discardDuplicates() {
  while (!pq.empty() && pq.top().first > dists[pq.top().second]) {
    pq.pop();
  }
}

IncompleteDijkstra::IncompleteDijkstra(const Graph* G, const std::vector<node>& sources,
                                       const std::unordered_set<node>* explored)
      : G(G), explored(explored) {
  if (!G) {
    throw invalid_argument("G is null");
  }

  for (node source : sources) {
    if (!explored || explored->find(source) == explored->end()) {
      dists[source] = 0.0;
      pq.emplace(0.0, source);
    }
  }
}

bool IncompleteDijkstra::hasNext() {
  discardDuplicates();
  return pq.size() > 0;
}

std::pair<node, edgeweight> IncompleteDijkstra::next() {
  if (!hasNext()) {
    throw std::invalid_argument("No next element");
  }

  // Extract nearest node
  discardDuplicates();
  node u;
  edgeweight dist_u;
  tie(dist_u, u) = pq.top();
  pq.pop();

  // Relax all of its edges
  G->forNeighborsOf(u, [&] (node v, edgeweight dist_uv) {
    if (explored && explored->find(v) != explored->end()) {
      return;
    }

    edgeweight new_dist_v = dist_u + dist_uv;
    auto it_dist_v = dists.find(v);
    if (it_dist_v == dists.end() || new_dist_v < it_dist_v->second) {
      dists[v] = new_dist_v;
      pq.emplace(new_dist_v, v);
    }
  });

  return {u, dist_u};
}

}

