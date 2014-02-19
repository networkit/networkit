/*
 * Diameter.h
 *
 *  Created on: 19.02.2014
 *      Author: Daniel Hoske, Christian Staudt
 */

#ifndef DIAMETER_H_
#define DIAMETER_H_

#include "../graph/Graph.h"

namespace NetworKit {

class Diameter {

public:

	/**
	 * TODO: documentation
	 */
	static std::pair<count, count> estimatedDiameterRange(const Graph& G, double error);
};

} /* namespace NetworKit */

#endif /* DIAMETER_H_ */
