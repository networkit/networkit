/*
 * DynAPSP.h
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe, Elisabetta Bergamini
 */

#ifndef DYNAPSP_H_
#define DYNAPSP_H_

#include "APSP.h"
#include "../dynamics/GraphEvent.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Dynamic APSP.
 */
class DynAPSP : public APSP {

public:

	/**
	 * Creates the object for @a G and source @a s.
	 *
	 * @param G The graph.
	 * @param s The source node.
	 * @param   storePredecessors   keep track of the lists of predecessors?
	 */
	DynAPSP(const Graph& G);

	/** Updates the distances after an event.*/
	void update(const std::vector<GraphEvent>& batch);

private:
	void dynamic_sssp(node root, std::vector<std::pair<node, edgeweight> > batch);
	std::vector<std::vector<edgeweight> > L;
	const edgeweight infDist = std::numeric_limits<edgeweight>::max();
};

} /* namespace NetworKit */

#endif /* DYNDIJKSTRA_H_ */
