/*
 * NodeNormalizedTriangleScore.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#ifndef NODENORMALIZEDTRIANGLESCORE_H
#define NODENORMALIZEDTRIANGLESCORE_H

#include "../edgescores/EdgeScore.h"

namespace NetworKit {

class NodeNormalizedTriangleScore : public EdgeScore<double> {

public:
	NodeNormalizedTriangleScore(const Graph& G, const std::vector<count>& triangles);
	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;
	virtual void run() override;

private:
	const std::vector<count>& triangles;

};

} // namespace NetworKit

#endif // NODENORMALIZEDTRIANGLESCORE_H
