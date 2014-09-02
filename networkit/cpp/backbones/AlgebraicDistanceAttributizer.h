/*
 * AlgebraicDistanceAttributizer.h
 *
 *  Created on: 01.09.2014
 *      Author: Christian Staudt
 */

#ifndef ALGEBRAICDISTANCEATTRIBUTIZER_H_
#define ALGEBRAICDISTANCEATTRIBUTIZER_H_

#include "BackboneCalculator.h"
#include "../distmeasures/AlgebraicDistance.h"

namespace NetworKit {

/**
 * EXPERIMENTAL
 */
class AlgebraicDistanceAttributizer : public AttributeGenerator<int, double> {

public:

	AlgebraicDistanceAttributizer(AlgebraicDistance& adObj);


	std::vector<double> getAttribute(const Graph& G, const std::vector<int>& attribute);

private:

	AlgebraicDistance& adObj;

};

}
/* namespace NetworKit */
#endif /* ALGEBRAICDISTANCEATTRIBUTIZER_H_ */
