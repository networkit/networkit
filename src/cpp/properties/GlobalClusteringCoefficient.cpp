/*
 * GlobalClusteringCoefficient.cpp
 *
 *  Created on: 12.11.2013
 */

#include "GlobalClusteringCoefficient.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

int uniformRandom(int max) {
  static int offset = 0;
  
  int currentMax = 1;
  int currentValue = 0;
  while(currentMax < max) {
    currentValue = currentValue * RAND_MAX + Aux::Random::integer();
    currentMax *= RAND_MAX;
  }
  int value = currentValue % max;
  return offset = (value + offset) % max;
}

unsigned int findIndex(const std::vector<int>& w, int v,
                       unsigned int lowerIdx, unsigned int upperIdx) {
  if(upperIdx - lowerIdx <= 1) {
    return lowerIdx;
  }
  int middleIdx = (upperIdx + lowerIdx) / 2;
  if(v >= w[middleIdx]) {
    return findIndex(w, v, middleIdx, upperIdx);
  } else {
    return findIndex(w, v, lowerIdx, middleIdx);
  }
}

double GlobalClusteringCoefficient::approximate(const Graph& G, int k) {
  const count n = G.numberOfNodes();
  
  std::vector<int> w(n + 1);
  int sum = 0;
  for(node i = 0; i < n; i++) {
    w[i] = sum;
    sum += (G.degree(i) * (G.degree(i) - 1)) / 2;
  }
  w[n] = sum;

  int l = 0;
  for(int i = 0; i < k; i++) {
    int r2 = uniformRandom(w[n]);
    node r = findIndex(w, r2, 0, n);
    node u = G.randomNeighbor(r);
    node w;
    do {
      w = G.randomNeighbor(r);
    } while (w == u);
    if(G.hasEdge(u, w)) {
      l++;
    }
  }

  return (double)l / (double)k;
}

} /* namespace NetworKit */
