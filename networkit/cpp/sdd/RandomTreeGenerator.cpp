/*
 * RandomTreeGenerator.cpp
 *
 *  Created on: 19.04.2014
 *      Author: dhoske
 */

#include <unordered_map>

#include "RandomTreeGenerator.h"
#include "../auxiliary/PrioQueue.h"

using namespace std;

namespace NetworKit {
namespace SDD {

RootedTree RandomTreeGenerator::prueferToTree(count n, node root, const vector<node>& pruefer, const vector<edgeweight>& weights) {
  Graph G(n, true);

  /* Determine node degrees. */
  vector<int> degrees(n, 1);
  for (node u: pruefer) {
    degrees[u]++;
  }
  Aux::PrioQueue<int, node> degrees_prio(degrees);

  /* Add edges from Pruefer sequence. */
  int edge_idx = 0;
  for (node u: pruefer) {
    /* Extract leaf u and decrease keys */
    node v;
    int deg;
    tie(deg, v) = degrees_prio.extractMin();
    assert(deg == 1);
    --degrees[v];
    --degrees[u];
    degrees_prio.decreaseKey(degrees[u], u); /* Prio Queue does not have getValue?? */

    /* Add new edge (u, v) */
    G.addEdge(u, v, weights[edge_idx++]);
  }

  /* Add last edge. */
  node u, v;
  int deg_u, deg_v;
  tie(deg_u, u) = degrees_prio.extractMin();
  tie(deg_v, v) = degrees_prio.extractMin();
  assert(deg_u == 1 && deg_v == 1);
  G.addEdge(u, v, weights[edge_idx++]);

  return move(RootedTree(G, root));
}

RootedTree RandomTreeGenerator::makeTree(count n, edgeweight lower, edgeweight upper) {
  assert(n > 2);

  /* Random Prüfer sequence */
  auto random_node = uniform_int_distribution<node>(0, n - 1);
  auto random_weight = uniform_real_distribution<edgeweight>(lower, upper);

  vector<node> pruefer(n - 2);
  vector<edgeweight> weights(n - 1);
  for (count i = 0; i < n - 1; ++i) {
    if (i < n - 2) {
      pruefer[i] = random_node(rand);
    }
    weights[i] = random_weight(rand);
  }
  return prueferToTree(n, random_node(rand), pruefer, weights);
}

}
}
