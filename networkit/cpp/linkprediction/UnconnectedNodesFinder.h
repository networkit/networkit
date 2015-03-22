/*
 * UnconnectedNodesFinder.h
 *
 *  Created on: 20.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef UNCONNECTEDNODESFINDER_H_
#define UNCONNECTEDNODESFINDER_H_

#include "../graph/Graph.h"

#include <utility>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * 
 */
class UnconnectedNodesFinder {
private:
  const Graph& G;

public:
  explicit UnconnectedNodesFinder(const Graph& G);

  std::vector<std::pair<node, node>> findAll(count k);

  std::vector<std::pair<node, node>> findFromNode(node u, count k);

};

} // namespace NetworKit

#endif /* UNCONNECTEDNODESFINDER_H_ */