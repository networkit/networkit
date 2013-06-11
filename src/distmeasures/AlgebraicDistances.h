/*
 * AlgebraicDistances.h
 *
 *  Created on: Jun 11, 2013
 *      Author: Henning
 */

#ifndef ALGEBRAICDISTANCES_H_
#define ALGEBRAICDISTANCES_H_

#include "../graph/Graph.h"
#include "../auxiliary/RandomProbability.h"
#include <vector>

namespace NetworKit {

class AlgebraicDistances {
protected:
	std::vector<std::vector<double> > loads; //<! loads[i]: vector of loads of length n for one system
	count numSystems;
	count numIters;
	const Graph& g;

	void randomInit();

public:
	AlgebraicDistances(const Graph& graph);
	virtual ~AlgebraicDistances();

	void preprocess(count numberSystems, count numberIterations, double omega);
	double algdist(node u, node v, count norm = 2) const;
};

} /* namespace NetworKit */

#endif /* ALGEBRAICDISTANCES_H_ */
