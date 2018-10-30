/*
 * Luby.h
 *
 *  Created on: 27.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef LUBY_H_
#define LUBY_H_

#include "IndependentSetFinder.h"


namespace NetworKit {

/**
 * @ingroup independentset
 * DEPRECATED: put into code archive as nobody seems to be using it
 * Luby's parallel independent set algorithm.
 */
class [[deprecated]]
Luby: public IndependentSetFinder {

public:

	// FIXME: check correctness of implementation
	std::vector<bool> run(const Graph& G) override;

	std::string toString() const override;
};

} /* namespace NetworKit */
#endif /* LUBY_H_ */
