/*
 * MultiscaleAttributizer.h
 *
 *  Created on: 20.06.2014
 *      Author: Gerd Lindner
 */

#ifndef MULTISCALEATTRIBUTIZER_H_
#define MULTISCALEATTRIBUTIZER_H_

#include "BackboneCalculator.h"
#include "gtest/gtest_prod.h"

namespace NetworKit {

/** 
 * Calculates the multiscale backbone for a given input graph.
 *
 * See "Extracting the multiscale backbone of complex weighted networks" by Serrano et al.
 */
class MultiscaleAttributizer : public AttributeGenerator {

public:

	/**
	 * Creates a new instance of the Multiscale Backbone calculator.
	 * @param alpha 		filter parameter
	 */
	MultiscaleAttributizer();

	EdgeAttribute getAttribute(const Graph& graph, const EdgeAttribute& attribute);

private:
	//Private helper functions
	double getProbability(count degree, edgeweight normalizedWeight);

	FRIEND_TEST(MultiscaleBackboneGTest, testSimpleMultiscaleBackbone);
};

}
/* namespace NetworKit */
#endif /* SIMMELIANBACKBONE_H_ */
