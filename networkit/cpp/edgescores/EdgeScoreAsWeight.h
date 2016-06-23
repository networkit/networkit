/*
 * EdgeScoreAsWeight.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef EDGESCOREASWEIGHT_H
#define EDGESCOREASWEIGHT_H

#include "../graph/Graph.h"

namespace NetworKit {

class EdgeScoreAsWeight {

public:
	EdgeScoreAsWeight(const Graph& G, const std::vector<double>& score, bool squared = false, edgeweight offset = 1, edgeweight factor = 1);
	Graph calculate();

private:
	const Graph& G;
	const std::vector<double>& score;
	bool squared;
	edgeweight offset;
	edgeweight factor;

};

} // namespace NetworKit

#endif // EDGESCOREASWEIGHT_H
