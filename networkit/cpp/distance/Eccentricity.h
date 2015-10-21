/*
 * Eccentricity.h
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#ifndef ECCENTRICITY_H_
#define ECCENTRICITY_H_

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup distance
 */
class Eccentricity {

public:

	/**
	 * TODO: documentation
	 */
	static std::pair<node, count> getValue(const Graph& G, node u);
};

} /* namespace NetworKit */

#endif /* ECCENTRICITY_H_ */
