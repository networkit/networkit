/*
 * EdgeScoreNormalizer.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef EDGESCORENORMALIZER_H
#define EDGESCORENORMALIZER_H

#include "../graph/Graph.h"
#include "EdgeScore.h"

namespace NetworKit {

template <typename A>
class EdgeScoreNormalizer : public EdgeScore<double> {

public:
	EdgeScoreNormalizer(const Graph &G, const std::vector<A> &score, bool invert = false, double lower = 0, double upper = 1.0);

	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;
	virtual void run() override;

private:
	const std::vector<A> &input;
	bool invert;
	double lower, upper;
};

}

#endif /* EDGESCORENORMALIZER_H */
