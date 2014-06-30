/*
 * IndependentSetFinder.h
 *
 *  Created on: 27.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef INDEPENDENTSETFINDER_H_
#define INDEPENDENTSETFINDER_H_

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup independentset
 * Abstract base class for independent set algorithms.
 */
class IndependentSetFinder {


public:

	/** Default destructor */
	virtual ~IndependentSetFinder() = default;

	/**
	 * Returns a boolean vector of length n where vec[v] is @c true iff v is in the independent sets.
	 * @param[in]	G	The graph.
	 * @return A boolean vector of length n.
	 */
	virtual std::vector<bool> run(const Graph& G) = 0;

	/**
	 * Get string representation of the algorithm.
	 * @return The string representation of the algorithm.
	 */
	virtual std::string toString() const;

	/**
	 * Checks whether a set is independent.
	 * @param set The set which is supposed to be independent.
	 * @param The graph.
	 * @return @c true iff @a set is independent.
	 */
	bool isIndependentSet(const std::vector<bool>& set, const Graph& G) const;

};

} /* namespace NetworKit */
#endif /* INDEPENDENTSETFINDER_H_ */
