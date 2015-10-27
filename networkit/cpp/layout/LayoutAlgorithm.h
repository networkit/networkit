/*
 * LayoutAlgorithm.h
 *
 *  Created on: May 20 2015
 *      Author: Christian Staudt
 */

#ifndef LAYOUTALGORITHM_H_
#define LAYOUTALGORITHM_H_

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup layout
 *
 * Base class for graph layout algorithms, i.e. algorithms that compute
 * a layout in form of 2D coordinates for nodes.
 */
class LayoutAlgorithm {


public:

	LayoutAlgorithm(const Graph& G) : G(G) {

	};

	virtual void run() = 0;

	virtual std::vector<std::pair<double, double> > getLayout() {
		return layout;
	};

protected:

	const Graph& G;
	std::vector<std::pair<double, double> > layout;

};

} /* namespace NetworKit */
#endif /* LAYOUTALGORITHM_H_ */
