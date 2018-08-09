/*
 * ErdosRenyiGenerator.cpp
 *
 *  Created on: 21.01.2014
 *      Author: Henning, Manuel Penschuck <networkit@manuel.jetzt>
 */

#include "ErdosRenyiGenerator.h"
#include "ErdosRenyiEnumerator.h"
#include "../graph/GraphBuilder.h"

namespace NetworKit {

ErdosRenyiGenerator::ErdosRenyiGenerator(count nNodes, double prob, bool directed) :
	nNodes{nNodes},	prob{prob}, directed{directed}
{}

Graph ErdosRenyiGenerator::generate() {
	GraphBuilder builder(nNodes, false, directed);

	{
		ErdosRenyiEnumeratorDefault impl(nNodes, prob, directed);
		impl.forEdgesParallel([&](int tid, node u, node v) {
			builder.addHalfEdge(u, v);
		});
	}

	return builder.toGraph(true, false);
}

} /* namespace NetworKit */
