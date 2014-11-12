/*
 * ChibaNishizekiTriangleCounter.h
 *
 *  Created on: 22.05.2014
 *      Author: Gerd Lindner
 */

#ifndef TRIANGLE_COUNTER_H_
#define TRIANGLE_COUNTER_H_

#include "../graph/Graph.h"
#include "AttributeGenerator.h"

namespace NetworKit {

/**
 * An implementation of the triangle counting algorithm by Ortmann et al.
 */
class TriangleCounter : public AttributeGenerator<int, int> {

public:

	std::vector<int> getAttribute(const Graph& graph, const std::vector<int>& attribute);
};

} /* namespace NetworKit */

#endif /* TRIANGLE_COUNTER_H_ */
