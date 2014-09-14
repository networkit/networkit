/*
 * ChungLu.h
 *
 *  Created on: Dec 23, 2013
 *      Author: Henning
 */

#ifndef CHUNGLU_H_
#define CHUNGLU_H_

#include "StaticDegreeSequenceGenerator.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

/** 
 * @ingroup generators
 * Given an arbitrary degree sequence, the Chung-Lu generative model
 * will produce a random graph with the same expected degree sequence. 
 * 
 * see Aiello, Chung, Lu: A Random Graph Model for Massive Graphs
 */
class ChungLuGenerator: public StaticDegreeSequenceGenerator {
protected:
	unsigned int sum_deg;
	count n;

public:
	ChungLuGenerator(const std::vector<unsigned int>& degreeSequence);

	/**
	 * Generates graph with expected degree sequence seq.
	 */
	virtual Graph generate();
};

} /* namespace NetworKit */
#endif /* CHUNGLU_H_ */
