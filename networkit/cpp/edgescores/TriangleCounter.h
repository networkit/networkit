/*
 * TriangleCounter.h
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner, Michael Hamann
 */

#ifndef TRIANGLE_COUNTER_H_
#define TRIANGLE_COUNTER_H_

#include "../graph/Graph.h"
#include "EdgeAttribute.h"

namespace NetworKit {

/**
 * An implementation of the triangle counting algorithm by Ortmann et al.
 */
class TriangleCounter : public EdgeAttribute<count> {

protected:
	const Graph& G;

public:

	TriangleCounter(const Graph& G);

	virtual std::vector<count> getAttribute() override;
};

} /* namespace NetworKit */

#endif /* TRIANGLE_COUNTER_H_ */
