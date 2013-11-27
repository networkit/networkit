/*
 * CoreDecomposition.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Henning
 */

#include "CoreDecomposition.h"

namespace NetworKit {

CoreDecomposition::CoreDecomposition() {
}

CoreDecomposition::~CoreDecomposition() {
}

struct NodeWithDegree {
  NodeWithDegree(node _v, count _degree) : v(_v), degree(_degree) {}
  node v;
  count degree;
};

std::vector<count> CoreDecomposition::run(const Graph& G) {
  const count n = G.numberOfNodes();
  std::vector<count> coreness = std::vector<count>(n, 0);

  //Knotengrad bestimmen
  std::vector<count> degree = std::vector<count>(n, 0);
  G.forEdges([&](node u, node v) {
    degree[u]++;
    degree[v]++;
  });

  // Knoten in entsprechende Liste stecken
  auto nodesByDegree = new std::list<NodeWithDegree>[n];
  std::vector<std::list<NodeWithDegree>::iterator> nodePointer;
  for(node v = 0; v < n; v++) {
    nodesByDegree[degree[v]].push_front(NodeWithDegree(v, degree[v]));
    nodePointer.push_back(nodesByDegree[degree[v]].begin());
  }

	index i = 1;
	Graph G2 = G;
  while (i - 1 < n) {
    while(!nodesByDegree[i - 1].empty()) {
      auto pv = nodesByDegree[i - 1].begin();
      coreness[pv->v] = i - 1;
      G2.forEdgesOf(pv->v, [&](node w, node u) {
        auto pu = nodePointer[u];
        if(pu->degree > i - 1) {
          auto &oldList = nodesByDegree[pu->degree];
          auto &newList = nodesByDegree[pu->degree - 1];
          newList.splice(newList.begin(), oldList, pu);
          pu->degree--;
        }
      });
      nodesByDegree[pv->degree].erase(pv);
    }
    i++;
	}
  delete[] nodesByDegree;

  return coreness;
}

} /* namespace NetworKit */
