/*
 * Matcher.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef MATCHER_H_
#define MATCHER_H_

#include "Matching.h"

namespace NetworKit {

/**
 * @ingroup matching
 * Abstract base class for matching algorithms.
 */
class Matcher {
protected:
	Graph& G;

public:
	/**
	 * Constructor.
	 * @param[in] G Graph for which matching is to be computed.
	 */
	Matcher(Graph& G);

	/** Default destructor */
	virtual ~Matcher() = default;

	/**
	 * Run the matching algorithm on the stored graph and return a matching.
	 * @return A matching of the stored graph.
	 */
	virtual Matching run() = 0;
};

} /* namespace NetworKit */
#endif /* MATCHER_H_ */
