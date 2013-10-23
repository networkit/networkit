/*
 * PLM2.cpp
 *
 *  Created on: 22.10.2013
 *      Author: cls
 */


#if 0

#include "PLM2.h"
#include "../auxiliary/Log.h"
#include "../coarsening/ClusterContracter.h"
#include "../coarsening/ClusteringProjector.h"

#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>



namespace NetworKit {


PLM2::PLM2(std::string par, double gamma) : parallelism(par), gamma(gamma), anyChange(false) {
}


PLM2::~PLM2() {
}


Clustering PLM2::pass(const Graph& G) {
	// init clustering to singletons
	count n = G.numberOfNodes();
	Clustering zeta(n);
	zeta.allToSingletons();

	// $\omega(E)$
	edgeweight total = G.totalEdgeWeight();

	// For each node we store a map that maps from cluster ID
	// to weight of edges to that cluster, this needs to be updated when a change occurs
	std::vector<tbb::concurrent_unordered_map<cluster, edgeweight> > incidenceWeight(n);

	G.parallelForNodes([&](node u) {
		G.forWeightedEdgesOf(u, [&](node u, node v, edgeweight w) {
					cluster C = zeta[v];
					if (u != v) {
						incidenceWeight[u][C] += w;
					}
				});
	});

	// modularity update formula for node moves
	// $$\Delta mod(u:\ C\to D)=\frac{\omega(u|D)-\omega(u|C\setminus v)}{\omega(E)}+\frac{2\cdot\vol(C\setminus u)\cdot\vol(u)-2\cdot\vol(D)\cdot\vol(u)}{4\cdot\omega(E)^{2}}$$

	// parts of formula follow
	tbb::concurrent_unordered_map<node, double> volNode(n);

	// calculate and store volume of each node
	G.parallelForNodes([&](node u) {
		volNode[u] += G.weightedDegree(u);
		volNode[u] += G.weight(u, u); // consider self-loop twice
		});

	tbb::concurrent_unordered_map<cluster, double> volCluster(n);
	// set volume for all singletons
	zeta.parallelForEntries([&](node u, cluster C) {
		volCluster[C] = volNode[u];
	});


	// $\vol(C \ {x})$ - volume of cluster C excluding node x
	auto volClusterMinusNode = [&](cluster C, node x) {
		double volC = 0.0;
		double volN = 0.0;
		volC = volCluster[C];
		if (zeta[x] == C) {
			volN = volNode[x];
			return volC - volN;
		} else {
			return volC;
		}
	};


	// $\omega(u | C \ u)$
	// TODO: function can be replaced by map lookup
	auto omegaCut = [&](node u, cluster C) {
		edgeweight w = 0.0;
		w = incidenceWeight[u][C];
		return w;
	};

	// difference in modularity when moving node u from cluster C to D
	auto deltaMod = [&](node u, cluster C, cluster D) {
			double volN = 0.0;
			volN = volNode[u];
			double delta = (omegaCut(u, D) - omegaCut(u, C)) / total + this->gamma * ((volClusterMinusNode(C, u) - volClusterMinusNode(D, u)) * volN) / (2 * total * total);
			return delta;
		};


	// begin pass
	int i = 0;
	bool change = false; // change in last iteration?
	do {
		i += 1;
		DEBUG("---> PLM2 pass: iteration # " << i << ", parallelism: " << this->parallelism);
		change = false; // is clustering stable?

		// try to improve modularity by moving a node to neighboring clusters
		auto moveNode = [&](node u) {
//			std::cout << u << " ";

			cluster best = none;
			cluster C = none;
			cluster D = none;
			double deltaBest = -0.5;
// #pragma omp atomic read
			C = zeta[u];

			G.forNeighborsOf(u, [&](node v) {
				D = zeta[v];
				if (D != C) { // consider only nodes in other clusters (and implicitly only nodes other than u)
					double delta = deltaMod(u, C, D);
					if (delta > deltaBest) {
						deltaBest = delta;
						best = D;
					}
				}
			});

			if (deltaBest > 0.0) { // if modularity improvement possible
				assert (best != C && best != none);// do not "move" to original cluster
				double volN = 0.0;

				// update weight of edges to incident clusters
				G.forWeightedNeighborsOf(u, [&](node v, edgeweight w) {
					incidenceWeight[v][best] += w;
					incidenceWeight[v][C] -= w;
				});
// #pragma omp atomic write
				zeta[u] = best; // move to best cluster
// #pragma omp atomic read
				volN = volNode[u];
#pragma omp atomic write
				change = true; // change to clustering has been made
#pragma omp atomic write
				this->anyChange = true; // indicate globally that clustering was modified

				// update the volume of the two clusters
				volCluster[C] -= volN;
				volCluster[best] += volN;
			}
		};

		// apply node movement according to parallelization strategy
		if (this->parallelism == "none") {
			G.forNodes(moveNode);
		} else if (this->parallelism == "simple") {
			G.parallelForNodes(moveNode);
		} else if (this->parallelism == "balanced") {
			G.balancedParallelForNodes(moveNode);
		} else {
			ERROR("unknown parallelization strategy: " << this->parallelism);
			exit(1);
		}

	} while (change);

	return zeta;

}



Clustering PLM2::run(Graph& G) {
	INFO("starting Louvain method");

	// sub-algorithms
	ClusterContracter contracter;
	ClusteringProjector projector;

	// hierarchies
	std::vector<std::pair<Graph, NodeMap<node> > > hierarchy; // hierarchy of graphs G^{i} and maps M^{i->i+1}
	std::vector<NodeMap<node> > maps; // hierarchy of maps M^{i->i+1}

	int h = -1; // finest hierarchy level
	bool done = false; //

	Graph* graph = &G;

	do {
		h += 1; // begin new hierarchy level
		DEBUG("Louvain hierarchy level " << h);
		// one local optimization pass
		DEBUG("starting Louvain pass");
		Clustering clustering = this->pass(*graph);
		if (this->anyChange) {
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
	} while (!done);

	DEBUG("starting projection");
	// project fine graph to result clustering
	Clustering result = projector.projectCoarseGraphToFinestClustering(*graph,
			G, maps);
	return result;
}


std::string PLM2::toString() const {
	std::stringstream strm;
	strm << "PLM2(" << "parallelism=" << this->parallelism << ")";
	return strm.str();
}

} /* namespace NetworKit */

#endif