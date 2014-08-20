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
 * Calculates the multiscale backbone attribute for a given graph. Each edge is
 * assigned the maximum filter value in [0,1] for which the edge will be contained
 * in the multiscale backbone.
 *
 * See "Extracting the multiscale backbone of complex weighted networks" by Serrano et al.
 */
class MultiscaleAttributizer : public AttributeGenerator<double, double> {

public:

	/**
	 * Creates a new instance of the Multiscale attributizer.
	 */
	MultiscaleAttributizer();

	~MultiscaleAttributizer() = default;

	std::vector<double> getAttribute(const Graph& graph, const std::vector<double>& attribute);

private:
	//Private helper functions
	double getProbability(count degree, edgeweight normalizedWeight);

	FRIEND_TEST(MultiscaleBackboneGTest, testSimpleMultiscaleBackbone);
};

}
/* namespace NetworKit */
#endif /* SIMMELIANBACKBONE_H_ */
