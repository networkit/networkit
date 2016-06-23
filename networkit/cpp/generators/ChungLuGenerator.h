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
 * see Chung, Lu: The average distances in random graphs with given expected degrees
 * and Chung, Lu: Connected Components in Random Graphs with Given Expected Degree Sequences.
 * Aiello, Chung, Lu: A Random Graph Model for Massive Graphs describes a different generative model
 * which is basically asymptotically equivalent but produces multi-graphs.
 */
class ChungLuGenerator: public StaticDegreeSequenceGenerator {
protected:
	count sum_deg;
	count n;

public:
	ChungLuGenerator(const std::vector<count>& degreeSequence);

	/**
	 * Generates graph with expected degree sequence seq.
	 */
	virtual Graph generate();
};

} /* namespace NetworKit */
#endif /* CHUNGLU_H_ */
