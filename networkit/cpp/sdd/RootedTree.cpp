/*
 * RootedTree.cpp
 *
 *  Created on: Apr 26, 2014
 *      Author: dhoske
 */

#include <algorithm>

#include "RootedTree.h"
#include "LCA.h"
#include "../properties/ConnectedComponents.h"

using namespace std;

namespace NetworKit {
namespace SDD {

bool isConnected(const Graph& G) {
  /* Connected and m = n - 1 */
  ConnectedComponents con(G);
  con.run();
  return con.numberOfComponents() == 1;
}

bool isTree(const Graph& T) {
  /* Connected and m = n - 1 */
  return isConnected(T) && T.numberOfEdges() == T.numberOfNodes() - 1;
}

bool RootedTree::operator==(const RootedTree& that) const {
  if (root != that.root || children.size() != that.children.size()) {
    return false;
  }

  /* Same tree even if the order of the children is different */
  for (index u = 0; u < children.size(); ++u) {
    if (children[u].size() != that.children[u].size()
        || !is_permutation(children[u].begin(), children[u].end(), that.children[u].begin())) {
      return false;
    }
  }
  return true;
}

RootedTree::RootedTree(const Graph &G, node root)
    : children(G.numberOfNodes()), parent(G.numberOfNodes()), root(root) {
  assert(0 <= root && root < G.numberOfNodes());
  assert(isTree(G));

  /* Traverse the tree from the root down and store adjacencies. */
  vector<bool> visited(G.numberOfNodes(), false);
  G.BFSfrom(root, [&] (node u, edgeweight) {
    visited[u] = true;
    G.forEdgesOf(u, [&] (node u, node v, edgeweight weight) {
      if (!visited[v]) {
        children[u].emplace_back(v);
        parent[v] = {u, weight};
      }
    });
  });
}

std::vector<edgeweight> computeReciprocalDepths(const RootedTree& T) {
  std::vector<edgeweight> out(T.numberOfNodes(), 0.0);
  T.forNodesPre([&] (node u) {
    if (u != T.getRoot()) {
      node v = T.getParent(u);
      out[u] = out[v] + 1./T.getWeight(v, u);
    }
  });
  return out;
}

std::vector<edgeweight> computeResistances(const RootedTree& T, const std::vector<Edge>& edges) {
  LCA lca(T);
  auto depths = computeReciprocalDepths(T);

  // Once you have the LCA structure and the depths, each query is easy to answer.
  std::vector<edgeweight> out(edges.size());
  for (index i = 0; i < edges.size(); ++i) {
    out[i] = depths[edges[i].u] + depths[edges[i].v] - 2*depths[lca.query(edges[i].u, edges[i].v)];
  }
  return out;
}

}
}
