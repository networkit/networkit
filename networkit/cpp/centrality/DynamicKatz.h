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
	count k;
	count maxdeg;
	bool groupOnly;
	
	// Nodes that have Katz score that only differ by this constant might appear swapped in the ranking.
	double rankTolerance;

public:
	/**
	 * Constructs a DynamicKatz object for the given Graph @a G. The damping factor is set to 1/(maxdeg + 1), where maxdeg is the maxmum degree in the graph.
	 *
	 * @param[in] G The graph.
	 * @param[in] k The number k for which we want to find the top-k nodes with highest Katz centrality
	 */
	DynamicKatz(const Graph& G, count k, bool groupOnly = false,
			double tolerance = 1e-9);

	virtual void run();

	/**
  * Updates the katz centralities after an edge insertion or deletion on the graph.
  *
  * @param event The edge insertions or deletion.
  */
	void update(const std::vector<GraphEvent> &events);
	
	void update(GraphEvent singleEvent) {
		std::vector<GraphEvent> events{singleEvent};
		update(events);
	}

	node top(count n = 0) {
		assert(activeRanking.size() > n);
		return activeRanking[n];
	}

	/**
	 * Returns the (upper) bound of the centrality of each node
	 */
	double bound(node v);
	
	/**
	 * Returns true if the bounds are sharp enough to rank two nodes against each other.
	 */
	bool areDistinguished(node u, node v);

private:
	/**
	 * Returns true if the bounds are sharp enough to rank two nodes against each other.
	 * Precondition: The first node appears higher in the current ranking the the second one.
	 */
	bool areCorrectlyRanked(node high, node low);

	/**
	 * Performs a single iteration of the algorithm.
	 */
	void doIteration();
	
	/**
	 * Returns true iff the ranking converged for the top-k.
	 */
	bool checkConvergence();

	std::vector<bool> isActive;
	std::vector<node> activeRanking;

	std::vector<double> baseData;
	std::vector<double> boundData;

public: // TODO: This is public because tests access it.
	std::vector<std::vector<count>> nPaths;
	count levelReached = 0;
};

} /* namespace NetworKit */
#endif /* DYNAMICKATZ_H_ */
