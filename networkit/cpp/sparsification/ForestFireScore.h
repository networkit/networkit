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
 * Experimental
 */
class ForestFireScore : public EdgeScore<double> {

public:

	ForestFireScore(const Graph& graph, double pf, double targetBurntRatio);
	
	virtual void run() override;

private:
	double pf;
	double targetBurntRatio;

};

}
/* namespace NetworKit */
#endif /* FORESTFIRESCORE_H_ */
