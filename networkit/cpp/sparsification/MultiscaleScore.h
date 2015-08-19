/*
 * MultiscaleScore.h
 *
 *  Created on: 20.06.2014
 *      Author: Gerd Lindner
 */

#ifndef MULTISCALESCORE_H_
#define MULTISCALESCORE_H_

#include "../edgescores/EdgeScore.h"

namespace NetworKit {

/**
 * Calculates the multiscale edge score for a given graph. Each edge is
 * assigned the maximum filter value in [0,1] for which the edge will be contained
 * in the multiscale backbone.
 *
 * See "Extracting the multiscale backbone of complex weighted networks" by Serrano et al.
 */
class MultiscaleScore : public EdgeScore<double> {

public:

	MultiscaleScore(const Graph& graph, const std::vector<double>& attribute);
	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;
	virtual void run() override;
	double getProbability(count degree, edgeweight normalizedWeight);

private:
	const std::vector<double>& attribute;
};

}
/* namespace NetworKit */

#endif /* MULTISCALESCORE_H_ */
