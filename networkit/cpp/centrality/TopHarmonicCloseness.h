/*
 * TopHarmonicCloseness.h
 *
 * Created on: 25.02.2018
 *		 Author: Eugenio Angriman
 */

#ifndef TOPHARMONICCLOSENESS_H_
#define TOPHARMONICCLOSENESS_H_

#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup centrality
 */
class TopHarmonicCloseness : public Algorithm {
public:
  /**
   * Finds the top k nodes with highest harmonic closeness centrality faster
   * than computing it for all nodes. The implementation is based on "Computing
   * Top-k Centrality Faster in Unweighted Graphs", Bergamini et al., ALENEX16.
   * The algorithms are based on two heuristics. We reccommend to use large_diam
   * = false for complex networks (or networks with small diameter) and
   * large_diam = true for street networks (or networks with large diameters).
   * Notice that the worst case running time of the algorithm is O(nm), where n
   * is the number of nodes and m is the number of edges. However, for most
   * real-world networks the empirical running time is O(m).
   *
   * @param G An unweighted graph.
   * @param k Number of nodes with highest harmonic closeness that have to be
   * found. For example, k = 10, the top 10 nodes with highest harmonic
   * closeness will be computed.
   * @param first_heu If true, the neighborhood-based lower bound is computed
   * and nodes are sorted according to it. If false, nodes are simply sorted by
   * degree.
   * @param sec_heu If true, BFSbound is re-computed at each iteration. If
   * false, BFScut is used
   */
  TopHarmonicCloseness(const Graph &G, count k = 1, bool first_heu = true,
                       bool sec_heu = true);

  /**
   * Computes top-k harmonic closeness on the graph passed in the constructor.
   */
  void run();

  /**
   * Returns a list with the scores of the k nodes with the highest harmonic
   * closeness.
   */
  std::vector<node> topkNodesList();

  /**
   * Returns a list with the scores of the k nodes with highest harmonic
   * closeness.
   */
  std::vector<edgeweight> topkScoresList();

  /**
   * Returns the ranking of the top k nodes with highest harmonic closeness and
   * their score.
   */
  std::vector<std::pair<node, edgeweight>> ranking();

protected:
  Graph G;
  count k;
  bool first_heu, sec_heu;
  std::vector<node> topk;
};
} // namespace NetworKit
#endif /* TOPHARMONICCLOSENESS_H_ */
