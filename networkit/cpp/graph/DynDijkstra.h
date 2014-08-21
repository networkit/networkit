/*
 * DynDijkstra.h
 *
 *  Created on: 21.07.2014
 *      Author: ebergamini
 */

#ifndef DYNDIJKSTRA_H_
#define DYNDIJKSTRA_H_

#include "DynSSSP.h"
#include "Dijkstra.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Dynamic Dijkstra.
 */
class DynDijkstra : public DynSSSP, public Dijkstra {

public:

	/**
	 * Creates the object for @a G and source @a s.
	 *
	 * @param G The graph.
	 * @param s The source node.
	 */
	DynDijkstra(const Graph& G, node s);


	/** Updates the distances after an event.*/
	void update(const std::vector<GraphEvent>& batch) override;

protected:
	std::vector<Color> color;

};

} /* namespace NetworKit */

#endif /* DYNDIJKSTRA_H_ */
