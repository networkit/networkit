/*
 * ChungLu.h
 *
 *  Created on: Dec 23, 2013
 *      Author: Henning
 */

#ifndef CHUNGLU_H_
#define CHUNGLU_H_

#include "StaticGraphGenerator.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

/** 
 * Given an arbitrary degree sequence, the Chung-Lu generative model
 * will produce a random graph with the same degree sequence. 
 * 
 * see Aiello, Chung, Lu: A Random Graph Model for Massive Graphs
 */
class ChungLuGenerator: public StaticGraphGenerator {
protected:
	std::vector<unsigned long long> seq;
	count sum_deg;
	count n;

public:
	ChungLuGenerator(const std::vector<unsigned long long>& degreeSequence);
	virtual ~ChungLuGenerator();

	/**
	 * Generates graph with expected degree sequence seq.
	 */
	virtual Graph generate();
};

} /* namespace NetworKit */
#endif /* CHUNGLU_H_ */
