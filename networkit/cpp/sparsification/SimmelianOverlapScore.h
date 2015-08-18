/*
 * SimmelianOverlapScore.h
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#ifndef SIMMELIANOVERLAPSCORE_H_
#define SIMMELIANOVERLAPSCORE_H_

#include "SimmelianScore.h"
#include <set>

namespace NetworKit {

/**
 * Calculates the Simmelian backbone (paramaetric variant) for a given input graph.
 */
class SimmelianOverlapScore : public SimmelianScore {

public:

	/**
	 * Creates a new instance of the parametric variant of the Simmelian Backbone calculator.
	 * @param maxRank 		the maximum rank that is considered for overlap calculation
	 */
	SimmelianOverlapScore(const Graph& graph, const std::vector<count>& triangles, count maxRank);

	virtual std::vector<double> scores() override;
	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;
	virtual void run() override;

private:
	count maxRank;
};

}
/* namespace NetworKit */
#endif /* SIMMELIANOVERLAPSCORE_H_ */
