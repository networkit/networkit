/*
 * DynDijkstra2.h
 *
 *  Created on: 24.07.2014
 *      Author: ebergamini
 */

#ifndef DYNDIJKSTRA_H_
#define DYNDIJKSTRA_H_

#include "DynSSSP.h"
#include "Dijkstra.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Dynamic Dijkstra based on FMN approach.
 */
class DynDijkstra2 : public DynSSSP, public Dijkstra {

public:

	/**
	 * Creates the object for @a G and source @a s.
	 *
	 * @param G The graph.
	 * @param s The source node.
	 */
	DynDijkstra2(const Graph& G, node s);

	/** Updates the distances after an event.*/
	void update(const std::vector<GraphEvent>& batch) override;

protected:
	std::vector<Color> color;
	std::vector<Aux::PrioQueue<edgeweight, node> > N;

};

} /* namespace NetworKit */

#endif /* DYNDIJKSTRA_H_ */
