/*
 * ChungLu.h
 *
 *  Created on: Dec 23, 2013
 *      Author: Henning
 */

#ifndef CHUNGLU_H_
#define CHUNGLU_H_

#include "StaticGraphGenerator.h"

namespace NetworKit {

class ChungLuGenerator: public StaticGraphGenerator {
protected:
	std::vector<count> seq;

public:
	ChungLuGenerator(const std::vector<count>& degreeSequence);
	virtual ~ChungLuGenerator();

	/**
	 * Generates graph with expected degree sequence seq.
	 */
	virtual Graph generate();
};

} /* namespace NetworKit */
#endif /* CHUNGLU_H_ */
