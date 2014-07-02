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

public:

	/** Default destructor */
	virtual ~Matcher() = default;

	/**
	 * Run the matching algorithm on Graph @a G and return a matching.
	 * @return A matching of @a G
	 */
	virtual Matching run(Graph& G) = 0;
};

} /* namespace NetworKit */
#endif /* MATCHER_H_ */
