/*
 * DynKatzCentrality.h
 *
 *  Created on: April 2018
 *      Author: Alexander van der Grinten
 *      based on code by Elisabetta Bergamini
 */

#ifndef DYNKATZCENTRALITY_H_
#define DYNKATZCENTRALITY_H_

#include "Centrality.h"
#include "../auxiliary/PrioQueue.h"
#include "../dynamics/GraphEvent.h"

namespace NetworKit {

/**
 * @ingroup centrality
 * Finds the top-k nodes with highest Katz centrality
 */
class DynKatzCentrality: public Centrality {
protected:
	double alpha; // damping
	count k;
	count maxdeg;
	bool groupOnly;
	
	// Nodes that have Katz score that only differ by this constant might appear swapped in the ranking.
	double rankTolerance;
public:
	bool useQueue = false;

public:
	/**
	 * Constructs a DynKatzCentrality object for the given Graph @a G. The damping factor is set to 1/(maxdeg + 1), where maxdeg is the maxmum degree in the graph.
	 *
	 * @param[in] G The graph.
	 * @param[in] k The number k for which we want to find the top-k nodes with highest Katz centrality
	 */
	DynKatzCentrality(const Graph& G, count k, bool groupOnly = false,
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
	 * Returns true if the bounds are sharp enough to rank two nodes against
	 * each other **within the tolerance**.
	 * Precondition: The first node appears higher in the current ranking the the second one.
	 */
	bool areSufficientlyRanked(node high, node low);

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
	Aux::PrioQueue<double, node> activeQueue;
	bool filledQueue = false;

	std::vector<double> baseData;
	std::vector<double> boundData;

public: // TODO: This is public because tests access it.
	std::vector<std::vector<count>> nPaths;
	count levelReached = 0;
};

} /* namespace NetworKit */
#endif /* DYNKATZCENTRALITY_H_ */
