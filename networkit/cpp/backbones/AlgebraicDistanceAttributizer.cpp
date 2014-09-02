/*
* AlgebraicDistanceAttributizer.h
*
*  Created on: 01.09.2014
*      Author: Christian Staudt
*/

#include "AlgebraicDistanceAttributizer.h"

namespace NetworKit {

AlgebraicDistanceAttributizer::AlgebraicDistanceAttributizer(AlgebraicDistance& adObj) : adObj(adObj) {}

std::vector<double> AlgebraicDistanceAttributizer::getAttribute(const Graph& G, const std::vector<int>& triangles) {
	// run algebraic distance calculation

	INFO("calculating algebraic distances");
	adObj.preprocess();

	INFO("attributizing edges");
	std::vector<double> ad(G.upperEdgeIdBound(), 0.0);

	G.parallelForEdges([&](node u, node v, edgeid uv){
		ad[uv] = adObj.distance(u, v);
	});

	return ad;
}

} /* namespace NetworKit */
