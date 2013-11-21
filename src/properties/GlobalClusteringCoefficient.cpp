/*
 * GlobalClusteringCoefficient.cpp
 *
 *  Created on: 12.11.2013
 */

#include "GlobalClusteringCoefficient.h"

namespace NetworKit {

double GlobalClusteringCoefficient::calculate(const Graph& G) {
  int numerator=0;
  int denominator=0;
  int numerator_vorher=0;
  
  G.forNodes([&](node u)
  {
    if(G.degree(u)>=2) {

      numerator_vorher=numerator;
    
      denominator+=(G.degree(u))*(G.degree(u)-1)/2;
      G.forEdgesOf(u,[&](node u, node w)
        {  
          
          G.forEdgesOf(w,[&](node w, node v)
          {
            if(v!=u)
            {
              if(G.hasEdge(v,u)==1)
                numerator++;// Jedes Dreieck wird von drei unterschiedlichen Knoten gezählt,                         // deshalb Faktor 3 weggelassen.
            }
          });
        });
    numerator=(numerator-numerator_vorher)/2+numerator_vorher;// Für jeden Knoten wird jedes Dreieck doppelt gezählt.
                    // Das muss rausgerechnet werden.
    }
  });
  double coefficient = (double) numerator / denominator;
  return coefficient;
}

int uniformRandom(int max) {
  static int offset = 0;
  
  int currentMax = 1;
  int currentValue = 0;
  while(currentMax < max) {
    currentValue = currentValue * RAND_MAX + rand();
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
  srand(time(NULL));
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
