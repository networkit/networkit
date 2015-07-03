/*
 * RandomEdgeAttributizer.h
 *
 *  Created on: 11.08.2014
 *      Author: Gerd Lindner
 */

#ifndef RANDOMEDGEATTRIBUTIZER_H_
#define RANDOMEDGEATTRIBUTIZER_H_

#include "../edgeattributes/EdgeAttribute.h"

namespace NetworKit {

/**
 * Generates a random edge attribute. Each edge is assigned a random value in [0,1].
 */
class RandomEdgeAttributizer : public EdgeAttribute<double> {

public:

	/**
	 * Creates a new instance of the Random edge attributizer.
	 */
	RandomEdgeAttributizer(const Graph& graph);

	virtual std::vector<double> getAttribute() override;

private:
	const Graph& graph;

};

}
/* namespace NetworKit */
#endif /* RANDOMEDGEATTRIBUTIZER_H_ */
