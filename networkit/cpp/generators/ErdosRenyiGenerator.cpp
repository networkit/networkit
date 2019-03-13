/*
 * ErdosRenyiGenerator.cpp
 *
 *  Created on: 21.01.2014
 *      Author: Henning, Manuel Penschuck <networkit@manuel.jetzt>
 */

#include "../../include/networkit/generators/ErdosRenyiGenerator.hpp"
#include "../../include/networkit/generators/ErdosRenyiEnumerator.hpp"
#include "../../include/networkit/graph/GraphBuilder.hpp"

namespace NetworKit {

ErdosRenyiGenerator::ErdosRenyiGenerator(count nNodes, double prob, bool directed, bool self_loops) :
	nNodes{nNodes},	prob{prob}, directed{directed}, self_loops{self_loops}
{}

Graph ErdosRenyiGenerator::generate() {
	GraphBuilder builder(nNodes, false, directed);

	{
		ErdosRenyiEnumeratorDefault impl(nNodes, prob, directed);
		impl.forEdgesParallel([&](int /*tid*/, node u, node v) {
			if (!self_loops && u == v) return;
			builder.addHalfEdge(u, v);
		});
	}

	return builder.toGraph(true, false);
}

} /* namespace NetworKit */
