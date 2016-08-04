/*
 * SSSP.h
 *
 *  Created on: 15.04.2014
 *      Author: cls
 */

#ifndef SSSP_H_
#define SSSP_H_

#include "APP.h"


namespace NetworKit {

/**
 * @ingroup graph
 * Abstract base class for single-source shortest path algorithms.
 */
class SSSP: public APP<edgeweight> {

public:

	/**
	 * Creates the SSSP class for @a G and source @a s.
	 *
	 * @param G The graph.
	 * @param s The source node.
	 */
	SSSP(const Graph& G, node s, bool storePaths=true, bool storeStack=false, node target = none) : APP<edgeweight>(G, s, storePaths, storeStack, target) {

    }

	virtual ~SSSP() = default;

	/** Computes the shortest paths from the source to all other nodes. */
	virtual void run() = 0;

protected:

    using APP<edgeweight>::G;
    using APP<edgeweight>::source;
    using APP<edgeweight>::target;
    using APP<edgeweight>::distances;
    using APP<edgeweight>::previous;
    using APP<edgeweight>::npaths;
    using APP<edgeweight>::stack;
    using APP<edgeweight>::storePaths;
    using APP<edgeweight>::storeStack;
    using APP<edgeweight>::edgeWeights;
};

} /* namespace NetworKit */

#endif /* SSSP_H_ */
