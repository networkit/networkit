/*
 * Louvain.cpp
 *
 *  Created on: 25.02.2013
 *      Author: cls
 */

#include "Louvain.h"

namespace EnsembleClustering {

Louvain::Louvain(std::string par) {
	this->parallelism = par;
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

	// For each node we store a map that maps from cluster ID
	// to weight of edges to that cluster, this needs to be updated when a change occurs
	std::vector<std::map<cluster, edgeweight> > incidenceWeight(G.numberOfNodes());
	G.forWeightedEdges([&](node u, node v, edgeweight w) { // FIXME: parallel would be better, but might cause problems
		cluster C = zeta[v];
		if (u != v) {
#pragma omp critical
			{
				incidenceWeight[u][C] += w;
			}
		}
	});

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
	auto omegaCut = [&](node u, cluster C) {
		edgeweight w = 0.0;
#pragma omp critical
		{
			w = incidenceWeight[u][C];
		}
		return w;

//		edgeweight sum = 0.0;
//		G.forWeightedEdgesOf(u, [&](node u, node v, edgeweight w){
//			if (zeta[v] == C) {
//				if (v != u){
//					sum += w;
//				}
//			}
//		});
//		return sum;
	};


	// difference in modularity when moving node u from cluster C to D
	auto deltaMod = [&](node u, cluster C, cluster D){
		double delta = (omegaCut(u, D) - omegaCut(u, C)) / total + ((volClusterMinusNode(C, u) - volClusterMinusNode(D, u)) * volNode[u]) / (2 * total * total);
		return delta;
	};


	Luby luby;	// independent set algorithm
	std::vector<bool> I;
	if (this->parallelism == "independent") {
		INFO("finding independent set");
		I = luby.run(G);
	}

	// FIXME: parallel unit tests lead to hanging execution (deadlocks or infinite loops?) sometimes

	// begin pass
	int i = 0;
	bool change; // change in last iteration?
	do {
		i += 1;
		DEBUG("---> Louvain pass: iteration # " << i);
		change = false; // is clustering stable?

		// try to improve modularity by moving a node to neighboring clusters
		auto moveNode = [&](node u){
			cluster C = zeta[u];
			TRACE("Processing node " << u << " of cluster " << C);
//			std::cout << ".";
			cluster best;
			double deltaBest = -0.5;
			G.forNeighborsOf(u, [&](node v){
				TRACE("Neighbor " << v << ", which is still in cluster " << zeta[v]);
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

				TRACE("Move vertex " << u << " to cluster " << best << ", deltaMod: " << deltaBest);

				// update weight of edges to incident clusters
				G.forWeightedNeighborsOf(u, [&](node v, edgeweight w) {
#pragma omp critical
					{
						incidenceWeight[v][zeta[u]] -= w;
						incidenceWeight[v][best] += w;
					}
				});

				zeta[u] = best; // move to best cluster
				change = true;	// change to clustering has been made
				this->anyChange = true;	// indicate globally that clustering was modified
				// update the volume of the two clusters
#pragma omp atomic update
				volCluster[C] -= volNode[u];
#pragma omp atomic update
				volCluster[best] += volNode[u];
			}
		};

		// apply node movement according to parallelization strategy
		if (this->parallelism == "none") {
			G.forNodes(moveNode);
		} else if (this->parallelism == "naive") {
			G.parallelForNodes(moveNode);
		} else if (this->parallelism == "naive-balanced") {
			G.balancedParallelForNodes(moveNode);
		} else if (this->parallelism == "independent") {
			// try to move only the nodes in independent set
			G.parallelForNodes([&](node u){
				if (I[u]) {
					moveNode(u);
				}
			});
		} else {
			ERROR("unknown parallelization strategy: " << this->parallelism);
			exit(1);
		}

//		std::cout << std::endl;
	} while (change && i < MAX_LOUVAIN_ITERATIONS);

	return zeta;
}

Clustering Louvain::run(Graph& G) {
	INFO("starting Louvain method");

	// sub-algorithms
	ClusterContracter contracter;
	ClusteringProjector projector;

	// hierarchies
	std::vector<std::pair<Graph, NodeMap<node> > > hierarchy; // hierarchy of graphs G^{i} and maps M^{i->i+1}
	std::vector<NodeMap<node> > maps; // hierarchy of maps M^{i->i+1}

	int h = -1; 	// finest hierarchy level
	bool done = false;	//

	Graph* graph = &G;

	do {
		h += 1; // begin new hierarchy level
		INFO("Louvain hierarchy level " << h);
		// one local optimization pass
		DEBUG("starting Louvain pass");
		Clustering clustering = this->pass(*graph);
		if (this->anyChange){
			// contract the graph according to clustering
			DEBUG("starting contraction");
			hierarchy.push_back(contracter.run(*graph, clustering));
			maps.push_back(hierarchy[h].second);
			graph = &hierarchy[h].first;

			done = false; // new hierarchy level, continue with loop
			this->anyChange = false; // reset change flag
		} else {
			done = true; // if clustering was not modified, do not contract and exit loop
		}
	} while (! done);

	DEBUG("starting projection");
	// project fine graph to result clustering
	Clustering result = projector.projectCoarseGraphToFinestClustering(*graph, G, maps);
	return result;
}

std::string Louvain::toString() const {
	std::stringstream strm;
	strm << "Louvain";
	return strm.str();
}

} /* namespace EnsembleClustering */
