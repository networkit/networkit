/*
 * ForestFireScore.h
 *
 *  Created on: 26.08.2014
 *      Author: Gerd Lindner
 */

#ifndef FORESTFIRESCORE_H_
#define FORESTFIRESCORE_H_

#include "../edgescores/EdgeScore.h"

namespace NetworKit {

/**
 * Based on the Forest Fire algorithm introduced by Leskovec et al.
 * The burn frequency of the edges is used as edge score.
 */
class ForestFireScore : public EdgeScore<double> {

public:

	ForestFireScore(const Graph& graph, double pf, double targetBurntRatio);
	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;
	virtual void run() override;

private:
	double pf;
	double targetBurntRatio;

};

}
/* namespace NetworKit */
#endif /* FORESTFIRESCORE_H_ */
