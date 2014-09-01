/*
 * TestAttributizer.cpp
 *
 *  Created on: 28.08.2014
 *      Author: Gerd Lindner
 */

#include "TestAttributizer.h"

namespace NetworKit {

TestAttributizer::TestAttributizer(count minDegree, double randomness) : minDegree(minDegree), randomness(randomness) {}

std::vector<double> TestAttributizer::getAttribute(const Graph& graph, const std::vector<int>& triangles) {
	std::vector<double> test(graph.upperEdgeIdBound(), 0.0);

	double maxt = (double) *std::max_element(std::begin(triangles), std::end(triangles));

	graph.forEdges([&](node u, node v, edgeid eid) {
		if (graph.degree(u) <= minDegree || graph.degree(v) <= minDegree) {
			test[eid] = 1.0;
		}

		double random = Aux::Random::probability();
		test[eid] = (randomness * random) + (1-randomness) * (triangles[eid] / maxt);
	});

	return test;
}

} /* namespace NetworKit */
