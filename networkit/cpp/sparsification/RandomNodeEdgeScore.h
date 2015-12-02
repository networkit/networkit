/*
 * RandomNodeEdgeScore.h
 *
 *  Created on: 20.11.2014
 *      Author: Michael Hamann
 */

#ifndef RANDOMNODEEDGESCORE_H
#define RANDOMNODEEDGESCORE_H

#include "../edgescores/EdgeScore.h"

namespace NetworKit {

class RandomNodeEdgeScore : public EdgeScore<double> {

public:
	RandomNodeEdgeScore(const Graph& graph, double rneRatio = 0.8);
	virtual void run() override;
	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;

private:
	double rneRatio;

};

} /* namespace NetworKit */

#endif // RANDOMNODEEDGESCORE_H
