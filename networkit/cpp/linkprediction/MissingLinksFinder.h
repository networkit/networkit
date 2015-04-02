/*
 * MissingLinksFinder.h
 *
 *  Created on: 20.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef MISSINGLINKSFINDER_H_
#define MISSINGLINKSFINDER_H_

#include "../graph/Graph.h"

#include <utility>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * 
 */
class MissingLinksFinder {
private:
  const Graph& G;

public:
  explicit MissingLinksFinder(const Graph& G);

  std::vector<std::pair<node, node>> findAll(count k);

  std::vector<std::pair<node, node>> findFromNode(node u, count k);

};

} // namespace NetworKit

#endif /* MISSINGLINKSFINDER_H_ */