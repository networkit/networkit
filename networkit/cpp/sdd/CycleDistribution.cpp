/*
 * CycleDistribution.cpp
 *
 *  Created on: Apr 26, 2014
 *      Author: dhoske
 */

#include "CycleDistribution.h"
#include "Flow.h"

using namespace std;

namespace NetworKit {
namespace SDD {

CycleDistribution::CycleDistribution(const Graph&, const RootedTree&,
                                     const std::vector<Edge>& off_tree_edges, const std::vector<edgeweight>&,
                                     SeedType seed)
    : engine(seed),
      off_tree_edges(off_tree_edges) {
  assert(off_tree_edges.size() > 0);
}

double CycleDistribution::getProb(const Edge& e) const {
  /* Standard assumption: edges are chosen uniformly at random. */
  if (find(off_tree_edges.begin(), off_tree_edges.end(), e) == off_tree_edges.end()) {
    return 0.0;
  } else {
    return 1.0;
  }
}

UniformCycleDistribution::UniformCycleDistribution(const Graph &G, const RootedTree &T,
                                                   const vector<Edge>& off_tree_edges, const std::vector<edgeweight>& stretches,
                                                   SeedType seed)
    : CycleDistribution(G, T, off_tree_edges, stretches, seed),
      random_edge(0, int(off_tree_edges.size()) - 1) {
}

StretchCycleDistribution::StretchCycleDistribution(const Graph &G, const RootedTree &T,
                                                   const vector<Edge>& off_tree_edges, const vector<edgeweight>& stretches,
                                                   SeedType seed)
    : CycleDistribution(G, T, off_tree_edges, stretches, seed),
      probs(stretches),
      random_edge(begin(probs), end(probs)) {
#if(DISCRETE_O1)
  max_prob = *max_element(begin(probs), end(probs));
#endif
}

/* Probability now weighed by stretch */
double StretchCycleDistribution::getProb(const Edge& e) const {
  auto id_it = find(off_tree_edges.begin(), off_tree_edges.end(), e);
  if (id_it == off_tree_edges.end()) {
    return 0.0;
  } else {
    return probs[distance(off_tree_edges.begin(), id_it)];
  }
}

}
}
