/*
 * ForestFireAttributizer.h
 *
 *  Created on: 26.08.2014
 *      Author: Gerd Lindner
 */

#ifndef FORESTFIREATTRIBUTIZER_H_
#define FORESTFIREATTRIBUTIZER_H_

#include "../edgescores/EdgeAttribute.h"

namespace NetworKit {

/**
 * Experimental
 */
class ForestFireAttributizer : public EdgeAttribute<double> {

public:

	ForestFireAttributizer(const Graph& graph, double pf, double targetBurntRatio);
	
	virtual std::vector<double> getAttribute() override;

private:
	const Graph& graph;
	double pf;
	double targetBurntRatio;

};

}
/* namespace NetworKit */
#endif /* FORESTFIREATTRIBUTIZER_H_ */
