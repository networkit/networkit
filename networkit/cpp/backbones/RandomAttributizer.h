/*
 * RandomAttributizer.h
 *
 *  Created on: 11.08.2014
 *      Author: Gerd Lindner
 */

#ifndef RANDOMATTRIBUTIZER_H_
#define RANDOMATTRIBUTIZER_H_

#include "BackboneCalculator.h"

namespace NetworKit {

/**
 * Generates a random edge attribute. Each edge is assigned a random value in [0,1].
 */
class RandomAttributizer : public AttributeGenerator<double, double> {

public:

	/**
	 * Creates a new instance of the Random edge attributizer.
	 */
	RandomAttributizer(double randomness);
	~RandomAttributizer() = default;

	std::vector<double> getAttribute(const Graph& graph, const std::vector<double>& attribute);

private:
	double randomness;

};

}
/* namespace NetworKit */
#endif /* RANDOMATTRIBUTIZER_H_ */
