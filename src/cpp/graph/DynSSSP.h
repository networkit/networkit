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

/**
 * @ingroup graph
 * Abstract base class for dynamic single-source shortest path algorithms.
 */
class DynSSSP : public SSSP {

public:

	/**
	 * Creates the DynSSSP class for @a G and source @a s.
	 *
	 * @param G The graph.
	 * @param s The source node.
	 */
	DynSSSP(const Graph& G, node s);

	/** Default destructor */
	virtual ~DynSSSP() = default;

	/** Computes the shortest paths from the source to all other nodes. */
	virtual void run() = 0;


    virtual void update(const std::vector<GraphEvent>& batch) = 0;

};

} /* namespace NetworKit */

#endif /* DYNSSSP_H_ */
