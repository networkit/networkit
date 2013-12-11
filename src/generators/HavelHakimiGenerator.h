/*
 * HavelHakimiGenerator.h
 *
 *  Created on: Dec 10, 2013
 *      Author: Henning
 */

#ifndef HAVELHAKIMIGENERATOR_H_
#define HAVELHAKIMIGENERATOR_H_

#include <vector>

#include "../graph/Graph.h"
#include "StaticGraphGenerator.h"

namespace NetworKit {

class HavelHakimiGenerator: public NetworKit::StaticGraphGenerator  {
private:
	std::vector<count> seq;
	bool realizable;

public:
	/**
	 * @param sequence Degree sequence to realize. Must be non-increasing.
	 */
	HavelHakimiGenerator(const std::vector<count>& sequence);
	virtual ~HavelHakimiGenerator();

	/**
	 * Erdšs-Gallai test if degree sequence seq is realizable.
	 */
	virtual bool isRealizable();

	/**
	 * Generates
	 */
	virtual Graph generate();
};

} /* namespace NetworKit */
#endif /* HAVELHAKIMIGENERATOR_H_ */
