/*
 * RandomEdgeScore.h
 *
 *  Created on: 11.08.2014
 *      Author: Gerd Lindner
 */

#ifndef RANDOMEDGESCORE_H_
#define RANDOMEDGESCORE_H_

#include "../edgescores/EdgeScore.h"

namespace NetworKit {

/**
 * Generates a random edge attribute. Each edge is assigned a random value in [0,1].
 */
class RandomEdgeScore : public EdgeScore<double> {

public:

	/**
	 * Creates a new instance of the Random edge score.
	 */
	RandomEdgeScore(const Graph& G);

	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;
	virtual void run() override;
};

}
/* namespace NetworKit */
#endif /* RANDOMEDGESCORE_H_ */
