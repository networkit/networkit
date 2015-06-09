/*
 * LCA.cpp
 *
 *  Created on: 19.07.2014
 *      Author: dhoske
 */

#include "LCA.h"

using namespace std;

namespace NetworKit {
namespace SDD {

LCA::LCA(const RootedTree& T) : n(T.numberOfNodes()), first_occ(n, index(-1)), depths(n) {
  // Transform from LCA into RMQ by Eulerian traversal
  std::vector<node> eulerian;
  T.forNodesEulerian([&] (node u, count depth) {
    if (first_occ[u] == index(-1)) {
      first_occ[u] = eulerian.size();
      depths[u] = depth;
    }
    eulerian.emplace_back(u);
  });

  // Build RMQ data structure
  rmq.emplace_back(eulerian);
  count size = eulerian.size();
  int idx = 1;
  // For all powers of two precompute RMQs
  for (count di = 1; di < size; di *= 2) {
    rmq.emplace_back(size);
    for (index i_left = 0; i_left < size; ++i_left) {
      index i_right = i_left + di;
      node min_left  = rmq[idx - 1][i_left];
      if (i_right < size) {
        node min_right = rmq[idx - 1][i_right];
        rmq[idx][i_left] = (depths[min_left] < depths[min_right]) ? min_left : min_right;
      } else {
        rmq[idx][i_left] = min_left;
      }
    }
    ++idx;
  }
}

node LCA::query(node u, node v) const {
  assert(u < n && v < n);

  // Subdivide query into two power-of-two queries
  index left  = min(first_occ[u], first_occ[v]);
  index right = max(first_occ[u], first_occ[v]);
  index level = static_cast<index>(floor(log2(static_cast<double>(right - left + 1))));

  // Answer the two power-of-two queries
  assert(level < rmq.size());
  node min_left  = rmq[level][left];
  node min_right = rmq[level][right - (1 << level) + 1];
  return (depths[min_left] < depths[min_right]) ? min_left : min_right;
}

}
}

