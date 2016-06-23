/*
 * TopCloseness.h
 *
 *  Created on: 03.10.2014
 *      Author: ebergamini, michele borassi
 */

#ifndef TOPCLOSENESS_H_
#define TOPCLOSENESS_H_
#include "../graph/Graph.h"
#include "../base/Algorithm.h"
#include "../auxiliary/PrioQueue.h"

namespace NetworKit {

/**
 * @ingroup centrality
 */
class TopCloseness : public Algorithm {
public:

  /**
	 * Finds the top k nodes with highest closeness centrality faster than computing it for all nodes, based on "Computing Top-k Closeness Centrality Faster in Unweighted Graphs", Bergamini et al., ALENEX16.
	 * The algorithms is based on two independent heuristics, described in the referenced paper. We recommend to use first_heu = true and second_heu = false for complex networks and first_heu = true and second_heu = true for street networks or networks with large diameters.
	 *
	 * @param G An unweighted graph.
	 * @param k Number of nodes with highest closeness that have to be found. For example, if k = 10, the top 10 nodes with highest closeness will be computed.
	 * @param first_heu If true, the neighborhood-based lower bound is computed and nodes are sorted according to it. If false, nodes are simply sorted by degree.
	 * @param sec_heu If true, the BFSbound is re-computed at each iteration. If false, BFScut is used.
	 * @
	 */
  TopCloseness(const Graph& G, count k = 1, bool first_heu = true, bool sec_heu = true);

  /**
	* Computes top-k closeness
	*
	*/
	void run();

  /**
	* Returns a list with the k nodes with highest closeness
	*
	*/
	std::vector<node> topkNodesList();

  /**
	* Returns a list with the scores of the k nodes with highest closeness
	*
	*/
	std::vector<edgeweight> topkScoresList();

protected:
	Graph G;
  count n;
	count k;
	bool first_heu, sec_heu;
	std::vector<node> topk;
	count visEdges = 0;
  count n_op = 0;
	std::vector<std::vector<node>> levels;
	std::vector<count> nodesPerLev;
  count nLevs = 0;
  std::vector<edgeweight> topkScores;
	std::vector<count> maxlevel;
	std::vector<count> maxlevelSize;
  std::vector<std::vector<count>> subtree;
  std::vector<double> farness;
  std::vector<count> reachL;
  std::vector<count> reachU;
  std::vector<count> component;

	void init();
    double BFScut(node v, double x, bool *visited, count *distances, node *pred, count *visEdges);
    void computelBound1(std::vector<double> &S);
    void BFSbound(node x, std::vector<double> &S, count *visEdges);
    void computeReachable();
    void computeReachableNodesUndir();
    void computeReachableNodesDir();
};


inline std::vector<node> TopCloseness::topkNodesList() {
	if (!hasRun) throw std::runtime_error("Call run method first");
	 return topk;
}

inline std::vector<edgeweight> TopCloseness::topkScoresList() {
	if (!hasRun) throw std::runtime_error("Call run method first");
	 return topkScores;
}

} /* namespace NetworKit */
#endif /* TOPCLOSENESS_H_ */
