/*
 * TestAttributizer.h
 *
 *  Created on: 28.08.2014
 *      Author: Gerd Lindner
 */

#ifndef TESTATTRIBUTIZER_H_
#define TESTATTRIBUTIZER_H_

#include "BackboneCalculator.h"

namespace NetworKit {

/**
 * EXPERIMENTAL
 */
class TestAttributizer : public AttributeGenerator<int, double> {

public:

	TestAttributizer(count minDegree, double randomness);
	~TestAttributizer() = default;

	std::vector<double> getAttribute(const Graph& graph, const std::vector<int>& attribute);

private:
	count minDegree;
	double randomness;

};

}
/* namespace NetworKit */
#endif /* TESTATTRIBUTIZER_H_ */
