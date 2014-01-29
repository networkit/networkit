/*
 * Louvain.cpp
 *
 *  Created on: 25.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu), Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#include "PLM.h"

#include "../coarsening/ClusterContracter.h"
#include "../coarsening/ClusteringProjector.h"
#include "omp.h"

namespace NetworKit {


PLM::PLM(std::string par, double gamma) : anyChange(false), parallelism(par), gamma(gamma) {

	this->VERSION = "1.0";
}

PLM::~PLM() {
	// TODO Auto-generated destructor stub
}

Clustering PLM::pass(Graph& G) {

	// FIXME: PLM cannot deal with deleted nodes

	// init clustering to singletons
	count n = G.numberOfNodes();
	Clustering zeta(n);
	zeta.allToSingletons();

	// $\omega(E)$
	edgeweight total = G.totalEdgeWeight();

	// For each node we store a map that maps from cluster ID
	// to weight of edges to that cluster, this needs to be updated when a change occurs
	std::vector<std::map<cluster, edgeweight> > incidenceWeight(
			G.numberOfNodes());
	G.parallelForNodes([&](node u) {
		G.forWeightedEdgesOf(u, [&](node u, node v, edgeweight w) {
			cluster C = zeta[v];
			if (u != v) {
				incidenceWeight[u][C] += w;
			}
		});
	});

#ifdef _OPENMP
	// create and init locks
	std::vector<omp_lock_t> mapLocks(n);
	G.parallelForNodes([&](node u) {
				omp_init_lock(&mapLocks[u]);
			});
#endif

	// modularity update formula for node moves
	// $$\Delta mod(u:\ C\to D)=\frac{\omega(u|D)-\omega(u|C\setminus v)}{\omega(E)}+\frac{2\cdot\vol(C\setminus u)\cdot\vol(u)-2\cdot\vol(D)\cdot\vol(u)}{4\cdot\omega(E)^{2}}$$

	// parts of formula follow
	std::vector<double> volNode(G.upperNodeIdBound(), 0.0);
	// calculate and store volume of each node
	G.parallelForNodes([&](node u) {
		volNode[u] += G.weightedDegree(u);
		volNode[u] += G.weight(u, u); // consider self-loop twice
		});

	std::vector<double> volCluster(G.upperNodeIdBound(), 0.0);
	// set volume for all singletons
	zeta.parallelForEntries([&](node u, cluster C) {
		volCluster[C] = volNode[u];
	});
	// end of initialization, set barrier for safety reasons

	// $\vol(C \ {x})$ - volume of cluster C excluding node x
	auto volClusterMinusNode = [&](cluster C, node x) {
		double volC = 0.0;
		double volN = 0.0;
// #pragma omp atomic read
		volC = volCluster[C];
		if (zeta[x] == C) {
			volN = volNode[x];
			return volC - volN;
		} else {
			return volC;
		}
	};

	// $\omega(u | C \ u)$
	auto omegaCut = [&](node u, cluster C) {
		edgeweight w = 0.0;
#ifdef _OPENMP
		omp_set_lock(&mapLocks[u]);
#endif
		#pragma omp atomic read
		w = incidenceWeight[u][C];
#ifdef _OPENMP
		omp_unset_lock(&mapLocks[u]);
#endif
		return w;
	};

	// difference in modularity when moving node u from cluster C to D
	auto deltaMod =
			[&](node u, cluster C, cluster D) {
			double volN = 0.0;
// #pragma omp atomic read
			volN = volNode[u];
			double delta = (omegaCut(u, D) - omegaCut(u, C)) / total + this->gamma * ((volClusterMinusNode(C, u) - volClusterMinusNode(D, u)) * volN) / (2 * total * total);
			return delta;
		};


	// begin pass
	int i = 0;
	bool change = false; // change in last iteration?
	do {
		i += 1;
		DEBUG("---> Louvain pass: iteration # " , i , ", parallelism: " , this->parallelism);
		change = false; // is clustering stable?

		// try to improve modularity by moving a node to neighboring clusters
		auto moveNode = [&](node u) {

			cluster best = none;
			cluster C = none;
			cluster D = none;
			double deltaBest = -0.5;

// #pragma omp atomic read
			C = zeta[u];

//			TRACE("Processing neighborhood of node " , u , ", which is in cluster " , C);
			G.forNeighborsOf(u, [&](node v) {
// #pragma omp atomic read
				D = zeta[v];
//				TRACE("Neighbor " , v , ", which is still in cluster " , zeta[v]);
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
#ifdef _OPENMP
					omp_set_lock(&mapLocks[v]);
#endif
					incidenceWeight[v][best] += w;
					incidenceWeight[v][C] -= w;
#ifdef _OPENMP
					omp_unset_lock(&mapLocks[v]);
#endif
				});
// #pragma omp atomic write
				zeta[u] = best; // move to best cluster
// #pragma omp atomic read
				volN = volNode[u];
// #pragma omp atomic write
				change = true; // change to clustering has been made
// #pragma omp atomic write
				this->anyChange = true; // indicate globally that clustering was modified

				// update the volume of the two clusters
#pragma omp atomic update
				volCluster[C] -= volN;
#pragma omp atomic update
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
			ERROR("unknown parallelization strategy: " , this->parallelism);
			exit(1);
		}

//		std::cout << std::endl;
	} while (change && i < MAX_LOUVAIN_ITERATIONS);

	// free lock mem
#ifdef _OPENMP
#pragma omp barrier
	DEBUG("about to destroy locks in Louvain pass");
	G.parallelForNodes([&](node u) {
				omp_destroy_lock(&mapLocks[u]);
			});
#endif

	return zeta;
}

Clustering PLM::run(Graph& G) {
	INFO("starting Louvain method");

	// sub-algorithms
	ClusterContracter contracter;
	ClusteringProjector projector;

	// hierarchies
	std::vector<std::pair<Graph, std::vector<node> > > hierarchy; // hierarchy of graphs G^{i} and maps M^{i->i+1}
	std::vector<std::vector<node> > maps; // hierarchy of maps M^{i->i+1}

	int h = -1; // finest hierarchy level
	bool done = false; //

	Graph* graph = &G;

	do {
		h += 1; // begin new hierarchy level
		DEBUG("Louvain hierarchy level " , h);
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

std::string PLM::toString() const {
	std::stringstream strm;
	strm << "PLM(" << this->parallelism << ")";
	return strm.str();
}

} /* namespace NetworKit */
