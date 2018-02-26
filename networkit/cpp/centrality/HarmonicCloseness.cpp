/*
 * HarmonicCloseness.cpp
 *
 * Created on: 24.02.2018
 * 		 Author: Eugenio Angriman
 */

#include <memory>

#include "../distance/BFS.h"
#include "../distance/Dijkstra.h"
#include "../distance/SSSP.h"
#include "HarmonicCloseness.h"

namespace NetworKit {

HarmonicCloseness::HarmonicCloseness(const Graph &G, bool normalized)
    : Centrality(G, normalized) {}

void HarmonicCloseness::run() {
  scoreData.assign(G.upperNodeIdBound(), 0.f);
  edgeweight infDist = std::numeric_limits<edgeweight>::max();

  G.parallelForNodes([&](node v) {
    std::unique_ptr<SSSP> sssp;
    if (G.isWeighted()) {
      sssp.reset(new Dijkstra(G, v, true, true));
    } else {
      sssp.reset(new BFS(G, v, true, true));
    }

    sssp->run();

    std::vector<edgeweight> distances = sssp->getDistances();

    double sum = 0;
    for (auto dist : distances) {
      if (dist != infDist && dist != 0) {
        sum += 1 / dist;
      }
    }

    scoreData[v] = sum;
  });
  if (normalized) {
    G.forNodes([&](node w) { scoreData[w] /= (G.numberOfNodes() - 1); });
  }

  hasRun = true;
}

double HarmonicCloseness::maximum() {
  return normalized ? (G.numberOfNodes() - 1) : 1.f;
}
} // namespace NetworKit
