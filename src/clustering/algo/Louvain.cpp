/*
 * Louvain.cpp
 *
 *  Created on: 25.02.2013
 *      Author: cls
 */

#include "Louvain.h"

namespace EnsembleClustering {

Louvain::Louvain() {
	// TODO Auto-generated constructor stub

}

Louvain::~Louvain() {
	// TODO Auto-generated destructor stub
}

Clustering Louvain::run(Graph& G) {

	Clustering zeta(G.numberOfNodes());
	zeta.allToSingletons();


	// parts of formula
	edgeweight total = G.totalEdgeWeight();

	// $\lamdbda(u)$
	auto lambda_1 = [&](node u){
		edgeweight sum = 0.0;
		G.forWeightedEdgesOf(u, [&](node u, node v, edgeweight w){
			sum += w;
		});
		sum += G.weight(u, u); // consider self-loop twice
		return sum;
	};

	// $\lambda(C)$
	auto lambda_2 = [&](cluster C){
		edgeweight sum = 0.0;
		G.forNodes([&](node u){
			if (zeta[u] == C) {
				sum += lambda_1(u);
			}
		});
		return sum;
	};


	// $\lambda(C \ x)$
	auto lambda_3 = [&](cluster C, node x){
		edgeweight sum = 0.0;
		G.forNodes([&](node u){
			if (zeta[u] == C) {
				if (u != x) {
					sum += lambda_1(u);
				}
			}
		});
		return sum;
	};

	auto omega_1 = [&](node u, cluster C){
		edgeweight sum = 0.0;
		G.forWeightedEdgesOf(u, [&](node u, node v, edgeweight w){
			if (zeta[v] == C) {
				sum += w;
			}
		});
		return sum;
	};

	// $\omega(u, C \ x)$
	auto omega_2 = [&](node u, cluster C){
		edgeweight sum = 0.0;
		G.forWeightedEdgesOf(u, [&](node u, node v, edgeweight w){
			if (zeta[v] == C) {
				if (v != u){
					sum += w;
				}
			}
		});
		return sum;
	};


	// difference in modularity when moving node u from cluster C to D
	auto deltaMod = [&](node u, cluster C, cluster D){
		double delta = (omega_1(u, D) - omega_2(u, C)) / total + (2 * lambda_3(C, u) * lambda_1(u) - 2 * lambda_2(D) * lambda_1(u)) / (4 * total * total);
		return delta;
	};


	bool stable = false;
	while (! stable){
		G.forNodes([&](node u){
			cluster C = zeta[u];
			cluster best;
			double deltaBest = -0.5;
			G.forNeighborsOf(u, [&](node v){
				cluster D = zeta[v];
				double delta = deltaMod(u, C, D);
				if (delta > deltaBest) {
					deltaBest = delta;
					best = D;
				}
				if (deltaBest > 0.0) {
					zeta.moveToCluster(best, u);
				}
			});
		});
	}


}

} /* namespace EnsembleClustering */
