/*
 *
 */

#include "AttributeAsWeight.h"

NetworKit::AttributeAsWeight::AttributeAsWeight(bool squared, edgeweight offset, edgeweight factor) : squared(squared), offset(offset), factor(factor) {
}

NetworKit::Graph NetworKit::AttributeAsWeight::calculate(NetworKit::Graph &g, const std::vector<double> &attribute) {
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
