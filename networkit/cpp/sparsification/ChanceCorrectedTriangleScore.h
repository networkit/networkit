/*
 * ChangeCorrectedTriangleScore.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#ifndef CHANCECORRECTEDTRIANGLESCORE_H
#define CHANCECORRECTEDTRIANGLESCORE_H

#include "../edgescores/EdgeScore.h"

namespace NetworKit {

class ChanceCorrectedTriangleScore : public EdgeScore<double> {

public:
	ChanceCorrectedTriangleScore(const Graph& graph, const std::vector<count>& triangles);
	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;
	virtual void run() override;

private:
	const std::vector<count>& triangles;

};

}
#endif // CHANCECORRECTEDTRIANGLESCORE_H
