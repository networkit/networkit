/*
 * DynAPSP.h
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe
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
	void update(const GraphEvent& event);

private:
	void dynamic_sssp(node source, node startbfs);
};

} /* namespace NetworKit */

#endif /* DYNDIJKSTRA_H_ */
