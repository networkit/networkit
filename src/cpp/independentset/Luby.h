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
 * Luby's parallel independent set algorithm.
 */
class Luby: public NetworKit::IndependentSetFinder {

public:

	// FIXME: check correctness of implementation
	std::vector<bool> run(const Graph& G) override;

	std::string toString() const override;
};

} /* namespace NetworKit */
#endif /* LUBY_H_ */
