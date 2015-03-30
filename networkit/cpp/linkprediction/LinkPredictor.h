/*
 * LinkPredictor.h
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef LINKPREDICTOR_H_
#define LINKPREDICTOR_H_

#include <memory>

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Abstract base class for link predictors.
 */
class LinkPredictor {
public:
  // Declare typedef in advance for later use
  typedef std::pair<std::pair<node, node>, double> node_dyad_score_pair;

private:
  /**
   * Implementation of the run method.
   * Doesn't have to check for a missing graph.
   *
   * @param u node in graph
   * @param v node in graph
   * @return a prediction-score indicating the likelihood of an
   * edge between the given nodes
   */
  virtual double runImpl(node u, node v) = 0;

protected:
  struct NodeDyadScoreComp {
    bool operator()(const node_dyad_score_pair& a, const node_dyad_score_pair& b) const {
      return (a.second > b.second) || (a.second == b.second && a.first < b.first);
    }
  } static ConcreteNodeDyadScoreComp;

  const Graph* G; //!< Graph to operate on

  bool validCache; //!< Indicates whether a possibly used cache is valid

public:
  explicit LinkPredictor();

  /**
   * Constructs a new LinkPredictor instance for the graph @a G.
   *
   * @param G The graph to work on
   */
  explicit LinkPredictor(const Graph& G);

  /**
   * Default destructor.
   */
  virtual ~LinkPredictor() = default;

  /**
   * Sets the graph to work on.
   * @param newGraph The graph to work on
   */
  void setGraph(const Graph& newGraph);

  /**
   * Predicts the likelihood of an edge between the given nodes.
   *
   * @param u node in graph
   * @param v node in graph
   * @return a prediction-score indicating the likelihood of an
   * edge between the given nodes
   */
  double run(node u, node v);

  virtual std::vector<node_dyad_score_pair> runOn(std::vector<std::pair<node, node>> nodePairs);

  virtual std::vector<node_dyad_score_pair> runOnParallel(std::vector<std::pair<node, node>> nodePairs);

  /**
   * Runs the link predictor on all node-pairs which are not connected
   * by a node in the given graph.
   *
   * @return a vector of dyad-score-pairs that is ordered descendingly by score and
   * on score equality ordered ascendingly by node-pairs.
   */
  virtual std::vector<node_dyad_score_pair> runAll();

  static std::vector<LinkPredictor::node_dyad_score_pair>& sortByScore(std::vector<node_dyad_score_pair>& predictions);
};

} // namespace NetworKit

#endif /* LINKPREDICTOR_H_ */