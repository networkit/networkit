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


Clustering Louvain::pass(Graph& G) {

	// init clustering to singletons
	Clustering zeta(G.numberOfNodes());
	zeta.allToSingletons();

	// $\omega(E)$
	edgeweight total = G.totalEdgeWeight();


	// modularity update formula for node moves
	// $$\Delta mod(u:\ C\to D)=\frac{\omega(u|D)-\omega(u|C\setminus v)}{\omega(E)}+\frac{2\cdot\vol(C\setminus u)\cdot\vol(u)-2\cdot\vol(D)\cdot\vol(u)}{4\cdot\omega(E)^{2}}$$

	// parts of formula follow
	NodeMap<double> volNode(G.numberOfNodes(), 0.0);
	// calculate and store volume of each node
	G.parallelForNodes([&](node u){
		volNode[u] += G.weightedDegree(u);
		volNode[u] += G.weight(u, u); // consider self-loop twice
	});

	IndexMap<cluster, double> volCluster(G.numberOfNodes(), 0.0);
	// set volume for all singletons
	zeta.parallelForEntries([&](node u, cluster C){
		volCluster[C] = volNode[u];
	});

	// $\vol(C \ {x})$ - volume of cluster C excluding node x
	auto volClusterMinusNode = [&](cluster C, node x){
		if (zeta[x] == C) {
			return volCluster[C] - volNode[x];
		} else {
			return volCluster[C];
		}

	};


	// $\omega(u | C \ u)$
	auto omegaCut = [&](node u, cluster C){
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
		double delta = (omegaCut(u, D) - omegaCut(u, C)) / total + (volClusterMinusNode(C, u) * volNode[u] - volClusterMinusNode(D, u) * volNode[u]) / (2 * total * total);
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
				if (zeta[v] != zeta[u]) { // consider only nodes in other clusters (and implicitly only nodes other than u)
					cluster D = zeta[v];
					double delta = deltaMod(u, C, D);
					if (delta > deltaBest) {
						deltaBest = delta;
						best = D;
					}
				}
			});
			if (deltaBest > 0.0) { // if modularity improvement possible
				assert (best != zeta[u]); // do not "move" to original cluster
				zeta[u] = best; // move to best cluster
				change = true;	// change to clustering has been made
				this->anyChange = true;	// indicate globally that clustering was modified
				// update the volume of the two clusters
				volCluster[C] -= volNode[u];
				volCluster[best] += volNode[u];
			}
		});
	} while (change);

	return zeta;
}

Clustering Louvain::run(Graph& G) {
	INFO("starting Louvain method");

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
		DEBUG("starting Louvain pass");
		clusterings.push_back(this->pass(graphs[h]));
		if (this->anyChange){
			// contract the graph according to clustering
			DEBUG("starting contraction");
			std::pair<Graph, NodeMap<node> > contraction = contracter.run(graphs[h], clusterings[h]);
			graphs.push_back(contraction.first);
			maps.push_back(contraction.second);

			done = false; // new hierarchy level, continue with loop
			this->anyChange = false; // reset change flag
		} else {
			done = true; // if clustering was not modified, do not contract and exit loop
		}
	} while (! done);

	DEBUG("starting projection");
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
