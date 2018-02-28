/*
 * DynTopHarmonicCloseness.h
 *
 *  Created on: 28.02.2018
 *      Author: nemes, Eugenio Angriman
 */

#ifndef DYNTOPHARMONICCLOSENESS_H_
#define DYNTOPHARMONICCLOSENESS_H_

#include "../base/Algorithm.h"
#include "../components/DynConnectedComponents.h"
#include "../components/DynWeaklyConnectedComponents.h"
#include "../graph/Graph.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 * Manages the list of the k nodes with the highest closeness centrality in a
 * dynamic graph.
 *
 * @ingroup centrality
 */
class DynTopHarmonicCloseness : public Algorithm {
public:
  /**
   * Constructs the DynTopHarmonicCloseness class. This class implements dynamic
   * algorithms for harmonic top-k closeness centrality. The implementation
   * is based on the static algorithms by Borassi et al. (complex networks)
   * and Bergamini et al. (large-diameter networks).
   *
   * @param G The graph.
   * @param k The number of top nodes.
   * @param useBFSbound Whether to use the algorithm for networks with large
   * diameter
   */
  DynTopHarmonicCloseness(const Graph &G, count k = 1,
                          bool useBFSbound = false);
  ~DynTopHarmonicCloseness();

  /**
   * Computes the k most central nodes on the initial graph.
   */
  void run();

  /**
   * Returns a list with the k most central nodes.
   * @return The list of the top-k nodes.
   */
  std::vector<node> topkNodesList();

  /**
   * Returns a list with the k highest closeness centralities.
   * @return The closeness centralities of the k most central nodes.
   */
  std::vector<edgeweight> topkScoresList();

  /**
   * Returns the ranking of the k most central nodes in the graph.
   * @return The ranking.
   */
  std::vector<std::pair<node, edgeweight>> ranking();

  void reset();

  /**
   * Updates the list of the k nodes with the highest closeness in G.
   *
   * @param event The edge modification event.
   */
  void update(const GraphEvent &event);

protected:
  /**
   * Runs a pruned BFS from the node @a u. The search is aborted once the
   * closeness centrality of @a u must be smaller than @a x.
   *
   * @param v The start node.
   * @param x The closeness centrality of the k-th most central node.
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
   * Computes the exact closeness centrality of the node @a x and stores upper
   * bounds for all other nodes in @a S2.
   *
   * @param x The start node.
   * @param S2 The upper bounds for all other nodes.
   * @param visEdges The number of visited edges.
   */
  void BFSbound(node x, std::vector<double> &S2, count *visEdges);

  /**
   * Handles an edge insertion.
   *
   * @param event The graph event.
   */
  void addEdge(const GraphEvent &event);
  /**
   * Handles an edge removal.
   *
   * @param event The graph event.
   */
  void removeEdge(const GraphEvent &event);

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

  /**
   * Recomputes the number of reachable nodes in an undirected graph after
   * inserting an edge between @a u and @a v.
   * @param u The first node.
   * @param v The second node.
   */
  void updateReachableNodesAfterInsertion(node u, node v);

  /**
   * Recomputes the number of reachable nodes in an undirected graph after
   * removing an edge between @a u and @a v.
   * @param u The first node.
   * @param v The second node.
   */
  void updateReachableNodesAfterDeletion(const GraphEvent &event);

  const Graph &G;
  count k = 1;
  bool useBFSbound;
  std::vector<node> topk;
  std::vector<edgeweight> topkScores;
  std::vector<edgeweight> allScores;
  std::vector<uint8_t> isExact;
  std::vector<uint8_t> isValid;
  std::vector<edgeweight> cutOff;
  std::vector<uint8_t> exactCutOff;
  DynConnectedComponents *comps;
  bool hasComps = false;
  Partition component;
  DynWeaklyConnectedComponents *wComps;
  bool hasWComps = false;
  std::vector<count> r;
  std::vector<count> rOld;
  std::vector<count> reachL;
};

inline std::vector<node> DynTopHarmonicCloseness::topkNodesList() {
  if (!hasRun) {
    throw std::runtime_error("Call run method first");
  }
  return topk;
}

inline std::vector<edgeweight> DynTopHarmonicCloseness::topkScoresList() {
  if (!hasRun) {
    throw std::runtime_error("Call run method first");
  }
  return topkScores;
}

} /* namespace NetworKit */
#endif /* DYNTOPCLOSENESSHARMONIC_H_ */
