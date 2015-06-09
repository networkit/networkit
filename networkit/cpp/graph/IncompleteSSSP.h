/*
 * IncompleteSSSP.h
 *
 *  Created on: 15.07.2014
 *      Author: dhoske
 */

#ifndef INCOMPLETESSSP_H_
#define INCOMPLETESSSP_H_

#include "Graph.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Abstract base class for single-source shortest path algorithms that return
 * the nodes in order of increasing distance from the source and do not
 * necessarily need to compute all distances.
 */
class IncompleteSSSP {

public:

  /**
   * Returns whether there is a next-nearest node
   * or all of the nodes reachable from the source
   * have already been processed.
   */
  virtual bool hasNext() = 0;

  /**
   * Returns the next-nearest node from the source and its
   * distance to the source. Should only be called if @a hasNext()
   * returns true.
   */
  virtual std::pair<node, edgeweight> next() = 0;
};

}

#endif
