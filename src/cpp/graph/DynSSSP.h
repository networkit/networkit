/*
 * DynSSSP.h
 *
 *  Created on: 17.07.2014
 *      Author: cls, ebergamini
 */

#ifndef DYNSSSP_H_
#define DYNSSSP_H_

#include <set>

#include "Graph.h"
#include "../dynamics/GraphEvent.h"
#include "SSSP.h"

namespace NetworKit {

  enum Color {WHITE, BLACK, GRAY};

/**
 * @ingroup graph
 * Interface for dynamic single-source shortest path algorithms.
 */
class DynSSSP {

public:

    virtual void update(const std::vector<GraphEvent>& batch) = 0;

};

} /* namespace NetworKit */

#endif /* DYNSSSP_H_ */
