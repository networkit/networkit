/*
 * ForestFireAttributizer.h
 *
 *  Created on: 26.08.2014
 *      Author: Gerd Lindner
 */

#ifndef FORESTFIREATTRIBUTIZER_H_
#define FORESTFIREATTRIBUTIZER_H_

#include "BackboneCalculator.h"

namespace NetworKit {

/**
 * Experimental
 */
class ForestFireAttributizer : public AttributeGenerator<int, double> {

public:

	ForestFireAttributizer(double pf, double targetBurntRatio);
	~ForestFireAttributizer() = default;
	std::vector<double> getAttribute(const Graph& graph, const std::vector<int>& attribute);

private:
	double pf;
	double targetBurntRatio;

};

}
/* namespace NetworKit */
#endif /* FORESTFIREATTRIBUTIZER_H_ */
