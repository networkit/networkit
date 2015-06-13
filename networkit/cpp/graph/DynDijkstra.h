/*
 * DynDijkstra.h
 *
 *  Created on: 21.07.2014
 *      Author: ebergamini
 */

#ifndef DYNDIJKSTRA_H_
#define DYNDIJKSTRA_H_

#include "DynSSSP.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Dynamic Dijkstra.
 */
class DynDijkstra : public DynSSSP {

public:

	/**
	 * Creates the object for @a G and source @a s.
	 *
	 * @param G The graph.
	 * @param s The source node.
	 * @param   storePredecessors   keep track of the lists of predecessors?
	 */
	DynDijkstra(const Graph& G, node s, bool storePredecessors = true);

	void run(node t = none) override;

	/** Updates the distances after an event.*/
	void update(const std::vector<GraphEvent>& batch) override;

protected:
	enum Color {WHITE, BLACK};
	std::vector<Color> color;

};


} /* namespace NetworKit */

#endif /* DYNDIJKSTRA_H_ */
