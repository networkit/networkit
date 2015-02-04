/*
 * RandomAttributizer.h
 *
 *  Created on: 11.08.2014
 *      Author: Gerd Lindner
 */

#ifndef RANDOMATTRIBUTIZER_H_
#define RANDOMATTRIBUTIZER_H_

#include "../edgeattributes/EdgeAttribute.h"

namespace NetworKit {

/**
 * Generates a random edge attribute. Each edge is assigned a random value in [0,1].
 */
class RandomAttributizer : public EdgeAttribute<double> {

public:

	/**
	 * Creates a new instance of the Random edge attributizer.
	 */
	RandomAttributizer(const Graph& graph);

	virtual std::vector<double> getAttribute() override;

private:
	const Graph& graph;

};

}
/* namespace NetworKit */
#endif /* RANDOMATTRIBUTIZER_H_ */
