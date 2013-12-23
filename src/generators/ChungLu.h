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

class ChungLu: public StaticGraphGenerator {
protected:
	std::vector<count> seq;

public:
	ChungLu(const std::vector<count>& degreeSequence);
	virtual ~ChungLu();

	/**
	 * Generates graph with expected degree sequence seq.
	 */
	virtual Graph generate();
};

} /* namespace NetworKit */
#endif /* CHUNGLU_H_ */
