/*
 * JarnikPrim.h
 *
 *  Created on: 13.05.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef JARNIKPRIM_H_
#define JARNIKPRIM_H_

#include "MSTFinder.h"
#include "../auxiliary/PrioQueue.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * The JarnikPrim class imlements the MST algorithm of Jarnik and Prim.
 */
class JarnikPrim : public MSTFinder {
public:
	/**
	 * Computes an MST for each connected component of @a G.
	 * @return A vector of MSTs, each being a vector of edges (node pairs).
	 */
	std::vector<std::vector<Edge>> run(const Graph &G);
};

} /* namespace NetworKit */

#endif /* JARNIKPRIM_H_ */
