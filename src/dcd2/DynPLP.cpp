/*
 * DynPLP.cpp
 *
 *  Created on: 03.01.2014
 *      Author: cls
 */

#include "DynPLP.h"

namespace NetworKit {

DynPLP::DynPLP(Graph& G, count theta) : DynCommunityDetector(G), updateThreshold(theta) {

}

void DynPLP::process(std::vector<GraphEvent>& stream) {
	auto isolate = [&](node u) {
		zeta[u] = zeta.addCluster();
	};

	for (GraphEvent ev : stream) {
		switch (ev.type) {
			case GraphEvent::NODE_ADDITION : {
				zeta.append(ev.u);
				zeta[ev.u] = zeta.addCluster();
				break;
			}
			case GraphEvent::NODE_REMOVAL : {
				zeta[ev.u] = none;
				break;
			}
			case GraphEvent::EDGE_ADDITION : {
				isolate(ev.u);
				isolate(ev.v);
				break;
			}
			case GraphEvent::EDGE_REMOVAL : {
				isolate(ev.u);
				isolate(ev.v);
				break;
			}
			case GraphEvent::EDGE_WEIGHT_UPDATE : {
				isolate(ev.u);
				isolate(ev.v);
				break;
			}
			case GraphEvent::TIME_STEP : {
				break;
			}
			default: {
				throw std::runtime_error("unknown event type");
			}
		}
	} // end event loop
}

Clustering DynPLP::retrieve() {
	typedef cluster label; // a label is the same as a cluster id
	count nUpdated; // number of nodes which have been updated in last iteration

	nIterations = 0; // number of iterations


	// propagate labels
	while (nUpdated > this->updateThreshold) {

		auto propagate = [&](node v) {
			if (G.degree(v) > 0) {
				std::map<label, double> labelWeights; // neighborLabelCounts maps label -> frequency in the neighbors

			}
		};

		G.balancedParallelForNodes(propagate);

	}


}

} /* namespace NetworKit */
