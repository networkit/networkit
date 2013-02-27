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


// TODO: optimize
Clustering Louvain::pass(Graph& G) {
	Clustering zeta(G.numberOfNodes());
	zeta.allToSingletons();

	// modularity update formula for node moves
	// $$\Delta mod(u:\ C\to D)=\frac{\omega(u|D)-\omega(u|C\setminus v)}{\omega(E)}+\frac{2\cdot\vol(C\setminus u)\cdot\vol(u)-2\cdot\vol(D)\cdot\vol(u)}{4\cdot\omega(E)^{2}}$$

	// parts of formula follow

	// $\omega(E)$
	edgeweight total = G.totalEdgeWeight();

	// $\lamdbda(u)$

	// TODO: call graph method
	auto vol_1 = [&](node u){
		edgeweight sum = 0.0;
		G.forWeightedEdgesOf(u, [&](node u, node v, edgeweight w){
			sum += w;
		});
		sum += G.weight(u, u); // consider self-loop twice
		return sum;
	};

	// $\vol(C)$
//	auto vol_2 = [&](cluster C){
//		edgeweight sum = 0.0;
//		G.forNodes([&](node u){
//			if (zeta[u] == C) {
//				sum += vol_1(u);
//			}
//		});
//		return sum;
//	};


	// $\vol(C \ x)$
	// TODO: call this volume
	auto vol_3 = [&](cluster C, node x){
		edgeweight sum = 0.0;
		G.forNodes([&](node u){ // FIXME: performance problem
			if (zeta[u] == C) {
				if (u != x) {
					sum += vol_1(u);
				}
			}
		});
		return sum;
	};

	// TODO: store and update volume for clusters - store vol_3 for every cluster - move one node - subtract wdge of v and add it to the other

	// $\omega(u | C)$
//	auto omega_1 = [&](node u, cluster C){
//		edgeweight sum = 0.0;
//		G.forWeightedEdgesOf(u, [&](node u, node v, edgeweight w){
//			if (zeta[v] == C) {
//				sum += w;
//			}
//		});
//		return sum;
//	};

	// $\omega(u | C \ u)$
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
		double delta = (omega_2(u, D) - omega_2(u, C)) / total + (vol_3(C, u) * vol_1(u) - vol_3(D, u) * vol_1(u)) / (2 * total * total);
		return delta;
	};

	int i = 0;
	bool change; // change in last iteration?
	do {
		i += 1;
		DEBUG("Louvain pass: iteration # " << i);
		change = false; // is clustering stable?
		G.forNodes([&](node u){
			cluster C = zeta[u];
			cluster best;
			double deltaBest = -0.5;
			G.forNeighborsOf(u, [&](node v){
				if (zeta[v] != zeta[u]) { // consider only nodes in other clusters
					cluster D = zeta[v];
					double delta = deltaMod(u, C, D);
					if (delta > deltaBest) {
						deltaBest = delta;
						best = D;
					}
				}
			});
			if ((deltaBest > 0.0)) { // if modularity improvement possible
				assert (best != zeta[u]); // do not "move" to original cluster
				zeta.moveToCluster(best, u);
				change = true;	// change to clustering has been made
				this->anyChange = true;	// indicate globally that clustering was modified
			}
		});
	} while (change);

	return zeta;
}

Clustering Louvain::run(Graph& G) {

	// sub-algorithms
	ClusterContracter contracter;
	ClusteringProjector projector;

	// hierarchies
	std::vector<Graph> graphs; // hierarchy of graphs G^{i}
	std::vector<Clustering> clusterings; // hierarchy of core clusterings \zeta^{i}
	std::vector<NodeMap<node> > maps; // hierarchy of maps M^{i->i+1}

	int h = -1; 	// finest hierarchy level
	bool done = false;	//

	graphs.push_back(G);

	do {
		h += 1; // begin new hierarchy level
		INFO("Louvain hierarchy level " << h);
		// one local optimization pass
		clusterings.push_back(this->pass(graphs[h]));
		if (this->anyChange){
			// contract the graph according to clustering
			std::pair<Graph, NodeMap<node> > contraction = contracter.run(graphs[h], clusterings[h]);
			graphs.push_back(contraction.first);
			maps.push_back(contraction.second);

			done = false; // new hierarchy level, continue with loop
			this->anyChange = false; // reset change flag
		} else {
			done = true; // if clustering was not modified, do not contract and exit loop
		}
	} while (! done);

	// project fine graph to result clustering
	Clustering result = projector.projectCoarseGraphToFinestClustering(graphs[h], graphs[0], maps);
	return result;
}

std::string Louvain::toString() const {
	std::stringstream strm;
	strm << "Louvain";
	return strm.str();
}

} /* namespace EnsembleClustering */
