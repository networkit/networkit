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
 * Abstract base class for graph coarsening/contraction algorithms.
 */
class GraphCoarsening {

public:

	GraphCoarsening();

	virtual ~GraphCoarsening();



};


} // namespace


#endif /* GRAPHCOARSENING_H_ */
