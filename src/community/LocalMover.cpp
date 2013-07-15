/*
 * LocalMover.cpp
 *
 *  Created on: 02.05.2013
 *      Author: cls
 */

#include "LocalMover.h"

namespace NetworKit {

LocalMover::LocalMover(QualityObjective& obj) {
	this->objective = &obj;
	// this->criterion = crit;
}

LocalMover::~LocalMover() {
	// TODO Auto-generated destructor stub
}

Clustering LocalMover::run(Graph& G) {
	Clustering zeta(G.numberOfNodes());

	this->zeta = &zeta;
	// TODO: initialize objective

	bool done; // TODO: define
	while (!done) {
		G.forNodes([&](node u) {
			// TODO: cluster C = objective->find(u);
			if (true) { // FIXME:
				// TODO: move(u, C);
			}
		});
	}

	return zeta;
}

void LocalMover::move(node u, cluster C) {
	(*zeta)[u] = C;
	// TODO: objective->onMove(u,C); // notify objective to update data structures
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


double NetworKit::LocalMover::Modularity::getValue(node v) {
}

double NetworKit::LocalMover::Coverage::getValue(node v) {
}
