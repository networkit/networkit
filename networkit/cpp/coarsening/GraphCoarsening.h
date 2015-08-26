/*
 * Contracter.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GRAPHCOARSENING_H_
#define GRAPHCOARSENING_H_

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup coarsening
 * Abstract base class for graph coarsening/contraction algorithms.
 */
class GraphCoarsening {

public:

	virtual ~GraphCoarsening() = default;

};


} // namespace


#endif /* GRAPHCOARSENING_H_ */
