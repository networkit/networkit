/*
 * LocalMover.cpp
 *
 *  Created on: 02.05.2013
 *      Author: cls
 */

#include "LocalMover.h"

namespace NetworKit {

LocalMover::LocalMover(Objective* obj, TerminationCriterion* crit) {
	this->objective = obj;
	this->criterion = crit;
}

LocalMover::~LocalMover() {
	// TODO Auto-generated destructor stub
}

Clustering LocalMover::run(Graph& G) {
	zeta = new Clustering(G.numberOfNodes());
	// TODO: initialize objective and termination criterion with graph and clustering

	while (! criterion->done()) {
		G.forNodes([&](node u) {
			cluster C = objective->find(u);
			if (true) { // FIXME:
				move(u, C);
			}
		});
	}

	return *zeta;
}

void LocalMover::move(node u, cluster C) {
	zeta[u] = C;
	objective->onMove(u,C); // notify objective to update data structures
}



//cluster NetworKit::LocalMover::DeltaModularity::find(node u) {
	// difference in modularity when moving node u from cluster C to D
//	auto deltaMod =
//			[&](node u, cluster C, cluster D) {
//			double volN = 0.0;
//#pragma omp atomic read
//			volN = volNode[u];
//			double delta = (omegaCut(u, D) - omegaCut(u, C)) / total + ((volClusterMinusNode(C, u) - volClusterMinusNode(D, u)) * volN) / (2 * total * total);
//			return delta;
//		};
//}


} /* namespace NetworKit */


