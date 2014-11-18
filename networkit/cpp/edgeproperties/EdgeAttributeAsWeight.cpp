/*
 * EdgeAttributeAsWeight.cpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#include "EdgeAttributeAsWeight.h"

namespace NetworKit {

AttributeAsWeight::AttributeAsWeight(bool squared, edgeweight offset, edgeweight factor) : squared(squared), offset(offset), factor(factor) {
}

Graph NetworKit::AttributeAsWeight::calculate(Graph &g, const std::vector<double> &attribute) {
	Graph result(g, true, false);

	if (squared) {
		g.parallelForEdges([&](node u, node v, edgeid eid) {
			result.setWeight(u, v, offset + factor * attribute[eid] * attribute[eid]);
		});
	} else {
		g.parallelForEdges([&](node u, node v, edgeid eid) {
			result.setWeight(u, v, offset + factor * attribute[eid]);
		});
	}

	return result;
}

} // namespace NetworKit
