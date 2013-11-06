/*
 * Luby.h
 *
 *  Created on: 27.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef LUBY_H_
#define LUBY_H_

#include "IndependentSetFinder.h"

#include "../auxiliary/RandomProbability.h"

namespace NetworKit {

class Luby: public NetworKit::IndependentSetFinder {

public:

	Luby() = default;

	virtual ~Luby() = default;

	std::vector<bool> run(const Graph& G) override;

	std::string toString() const override;
};

} /* namespace NetworKit */
#endif /* LUBY_H_ */
