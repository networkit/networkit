/*
 * DynamicKatz.h
 *
 *  Created on: 09.01.2015
 *      Author: Henning
 */

#ifndef DYNAMICKATZ_H_
#define DYNAMICKATZ_H_

#include "Centrality.h"
#include "../dynamics/GraphEvent.h"

namespace NetworKit {

/**
 * @ingroup centrality
 * Finds the top-k nodes with highest Katz centrality
 */
class DynamicKatz: public Centrality {
protected:
	double alpha; // damping
	count k = 100;
	count maxdeg;
	bool groupOnly;

public:
	/**
	 * Constructs a DynamicKatz object for the given Graph @a G. The damping factor is set to 1/(maxdeg + 1), where maxdeg is the maxmum degree in the graph.
	 *
	 * @param[in] G The graph.
	 * @param[in] k The number k for which we want to find the top-k nodes with highest Katz centrality
	 */
	DynamicKatz(const Graph& G, count k, bool groupOnly = true);

	virtual void run();

	/**
  * Updates the katz centralities after an edge insertion or deletion on the graph.
  *
  * @param event The edge insertions or deletion.
  */
	void update(GraphEvent event);

	// TODO make it private!!!
	std::vector<std::vector<count>> nPaths;
	count levelReached = 0;
};

} /* namespace NetworKit */
#endif /* DYNAMICKATZ_H_ */
