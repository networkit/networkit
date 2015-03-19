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
  struct NodeDyadScoreComp {
    bool operator()(const node_dyad_score_pair& a, const node_dyad_score_pair& b) const {
      return (a.second > b.second) || (a.second == b.second && a.first < b.first);
    }
  };

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

  /**
   * Implementation of the runAll method.
   * Doesn't have to check for a missing graph.
   *
   * @param limit Limit for the number of dyad-score-pairs to return.
   * If set to 0 all pairs will get returned.
   * @return a vector of dyad-score-pairs that is ordered descendingly by score and
   * on score equality ordered ascendingly by node-pairs.
   */
  virtual std::vector<node_dyad_score_pair> runAllImpl(count limit);

protected:
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

  /**
   * Runs the link predictor on all node-pairs which are not connected
   * by a node in the given graph.
   *
   * @param limit Limit for the number of dyad-score-pairs to return.
   * If set to 0 all pairs will get returned.
   * @return a vector of dyad-score-pairs that is ordered descendingly by score and
   * on score equality ordered ascendingly by node-pairs.
   */
  virtual std::vector<node_dyad_score_pair> runAll(count limit = 0);
};

} // namespace NetworKit

#endif /* LINKPREDICTOR_H_ */