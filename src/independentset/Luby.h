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

	Luby();

	virtual ~Luby();

	virtual std::vector<bool> run(const Graph& G);

	virtual std::string toString() const;
};

} /* namespace NetworKit */
#endif /* LUBY_H_ */
