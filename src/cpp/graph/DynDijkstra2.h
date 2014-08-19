/*
 * DynDijkstra2.h
 *
 *  Created on: 24.07.2014
 *      Author: ebergamini
 */

#ifndef DYNDIJKSTRA2_H_
#define DYNDIJKSTRA2_H_

#include "DynSSSP.h"
#include "Dijkstra.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Dynamic Dijkstra based on FMN approach.
 */
class DynDijkstra2 : public DynSSSP {

public:

	/**
	* Creates the object for @a G and source @a s.
	*
	* @param G The graph.
	* @param s The source node.
	* @param   storePredecessors   keep track of the lists of predecessors?
	*/
	DynDijkstra2(const Graph& G, node s, bool storePredecessors = true);

	void run(node t=none) override;

	/** Updates the distances after an event.*/
	void update(const std::vector<GraphEvent>& batch) override;

protected:
	enum Color {WHITE, BLACK};
	std::vector<Color> color;
	std::vector<Aux::PrioQueue<double, node> > N;
};

} /* namespace NetworKit */

#endif /* DYNDIJKSTRA2_H_ */
