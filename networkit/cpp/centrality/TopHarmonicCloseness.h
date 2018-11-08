/*
 * TopHarmonicCloseness.h
 *
 * Created on: 25.02.2018
 *		Authors: nemes, Eugenio Angriman
 */

#ifndef TOPHARMONICCLOSENESS_H_
#define TOPHARMONICCLOSENESS_H_

#include "../base/Algorithm.h"
#include "../components/ConnectedComponents.h"
#include "../graph/Graph.h"
#include "../structures/Partition.h"

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
   * The algorithms are based on two heuristics. We recommend to use
   * useBFSbound = false for complex networks (or networks with small diameter)
   * and useBFSbound = true for street networks (or networks with large
   * diameters). Notice that the worst case running time of the algorithm is
   * O(nm), where n is the number of nodes and m is the number of edges.
   * However, for most real-world networks the empirical running time is O(m).
   *
   * @param G An unweighted graph.
   * @param k Number of nodes with highest harmonic closeness that have to be
   * found. For example, k = 10, the top 10 nodes with highest harmonic
   * closeness will be computed.
   * @param useBFSbound If true, BFSbound is re-computed at each iteration. If
   * false, BFScut is used
   */
  TopHarmonicCloseness(const Graph &G, count k = 1, bool useBFSbound = false);

  /**
   * Computes top-k harmonic closeness on the graph passed in the constructor.
   */
  void run();

  /**
   * Returns a list with the k nodes with highest harmonic closeness.
   * WARNING: closeness centrality of some nodes below the top-k could be equal
   * to the k-th closeness, we call them trail. Set the parameter includeTrail
   * to true to also include those nodes but consider that the resulting vector
   * could be longer than k.
   *
   * @param includeTrail Whether or not to include trail nodes.
   * @return The list of the top-k nodes.
   */
  std::vector<node> topkNodesList(bool includeTrail = false);

  /**
   * Returns a list with the scores of the k nodes with highest harmonic
   * closeness.
   * WARNING: closeness centrality of some nodes below the top-k could
   * be equal to the k-th closeness, we call them trail. Set the parameter
   * includeTrail to true to also include those centrality values but consider
   * that the resulting vector could be longer than k.
   *
   * @param includeTrail Whether or not to include trail centrality value.
   * @return The closeness centralities of the k most central nodes.
   */
  std::vector<edgeweight> topkScoresList(bool includeTrail = false);

protected:
  /**
   * Runs a pruned BFS from the node @a u. The search is aborted once the
   * closeness centrality of @a u must be smaller than @a x.
   *
   * @param v The start node.
   * @param x The harmonic closeness centrality of the k-th most central node.
   * @param n The number of nodes in the graph.
   * @param r The number of reachable nodes (or an upper bound in directed
   * graphs) for each node.
   * @param visited The vector containing the visited status.
   * @param distances The vector containing the distances.
   * @param pred The vector containing the predecessor of each node.
   * @param visitedEdges The number of visited edges.
   * @return
   */
  std::pair<edgeweight, bool> BFScut(node v, edgeweight x, count n, count r,
                                     std::vector<uint8_t> &visited,
                                     std::vector<count> &distances,
                                     std::vector<node> &pred,
                                     count &visitedEdges);

  /**
   * Computes the exact harmonic closeness centrality of the node @a x and
   * stores upper bounds for all other nodes in @a S2.
   *
   * @param x The start node.
   * @param S2 The upper bounds for all other nodes.
   * @param visEdges The number of visited edges.
   */
  void BFSbound(node x, std::vector<double> &S2, count *visEdges);

  /**
   * Computes the number of reachable nodes or an upper bound in directed
   * graphs.
   */
  void computeReachableNodes();

  /**
   * Computes an upper bound for the number of reachable nodes for each node
   * in a directed graph.
   */
  void computeReachableNodesDirected();

  /**
   * Computes the number of reachable nodes in an undirected graph.
   */
  void computeReachableNodesUndirected();

  void init();

  Graph G;
  count k;
  count n;
  count trail = 0;
  double minCloseness = std::numeric_limits<double>::max();
  count nMinCloseness = 0;
  bool useBFSbound;
  std::vector<node> topk;
  std::vector<edgeweight> topkScores;
  std::vector<edgeweight> allScores;
  std::vector<uint8_t> isExact;
  std::vector<uint8_t> isValid;
  std::vector<edgeweight> cutOff;
  std::vector<uint8_t> exactCutOff;
  Partition component;
  std::vector<count> r;
  std::vector<count> rOld;
};

inline std::vector<node>
TopHarmonicCloseness::topkNodesList(bool includeTrail) {
  assureFinished();
  if (!includeTrail) {
    auto begin = topk.begin();
    std::vector<node> topkNoTrail(begin, begin + k);
    return topkNoTrail;
  }

  return topk;
}

inline std::vector<edgeweight>
TopHarmonicCloseness::topkScoresList(bool includeTrail) {
  assureFinished();
  if (!includeTrail) {
    auto begin = topkScores.begin();
    std::vector<double> topkScoresNoTrail(begin, begin + k);
    return topkScoresNoTrail;
  }
  return topkScores;
}

} // namespace NetworKit
#endif /* TOPHARMONICCLOSENESS_H_ */
