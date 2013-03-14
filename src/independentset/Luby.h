/*
 * Luby.h
 *
 *  Created on: 27.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef LUBY_H_
#define LUBY_H_

#include "IndependentSetFinder.h"

#include "../aux/RandomProbability.h"

namespace EnsembleClustering {

class Luby: public EnsembleClustering::IndependentSetFinder {

public:

	Luby();

	virtual ~Luby();

	virtual std::vector<bool> run(const Graph& G);

	virtual std::string toString() const;
};

} /* namespace EnsembleClustering */
#endif /* LUBY_H_ */
