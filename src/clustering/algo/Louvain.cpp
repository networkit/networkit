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
	const int NUM_THREADS = omp_get_num_threads();

#pragma omp barrier
	// create and init locks
	std::vector<omp_lock_t> mapLocks(n);
	G.parallelForNodes([&](node u) {
				omp_init_lock(&mapLocks[u]);
			});
#pragma omp barrier
#endif

	// modularity update formula for node moves
	// $$\Delta mod(u:\ C\to D)=\frac{\omega(u|D)-\omega(u|C\setminus v)}{\omega(E)}+\frac{2\cdot\vol(C\setminus u)\cdot\vol(u)-2\cdot\vol(D)\cdot\vol(u)}{4\cdot\omega(E)^{2}}$$

	// parts of formula follow
	NodeMap<double> volNode(G.numberOfNodes(), 0.0);
	// calculate and store volume of each node
	G.parallelForNodes([&](node u) {
		volNode[u] += G.weightedDegree(u);
		volNode[u] += G.weight(u, u); // consider self-loop twice
		});

	IndexMap<cluster, double> volCluster(G.numberOfNodes(), 0.0);
	// set volume for all singletons
	zeta.parallelForEntries([&](node u, cluster C) {
		volCluster[C] = volNode[u];
	});
	// end of initialization, set barrier for safety reasons
#pragma omp barrier

	// $\vol(C \ {x})$ - volume of cluster C excluding node x
	auto volClusterMinusNode = [&](cluster C, node x) {
		double volC = 0.0;
		double volN = 0.0;
//#ifdef _OPENMP
//		DEBUG("Try to acquire lock for pos " << C);
//		omp_set_lock(&mapLocks[C]);
//#endif
#pragma omp atomic read
		volC = volCluster[C];
//#ifdef _OPENMP
//		omp_unset_lock(&mapLocks[C]);
//		DEBUG("Released lock for pos " << C);
//#endif
		if (zeta[x] == C) {
//#ifdef _OPENMP
//			DEBUG("Try to acquire lock for pos " << x);
//			omp_set_lock(&mapLocks[x]);
//#endif
#pragma omp atomic read
			volN = volNode[x];
//#ifdef _OPENMP
//			omp_unset_lock(&mapLocks[x]);
//			DEBUG("Released lock for pos " << x);
//#endif
			return volC - volN;
		} else {
			return volC;
		}
	};

	// $\omega(u | C \ u)$
	auto omegaCut = [&](node u, cluster C) {
		edgeweight w = 0.0;
#ifdef _OPENMP
		DEBUG("Try to acquire lock for pos " << u);
			omp_set_lock(&mapLocks[u]);
#endif
			w = incidenceWeight[u][C];
#ifdef _OPENMP
			omp_unset_lock(&mapLocks[u]);
			DEBUG("Released lock for pos " << u);
#endif
			return w;
		};

	// difference in modularity when moving node u from cluster C to D
	auto deltaMod =
			[&](node u, cluster C, cluster D) {
//#ifdef _OPENMP
//		DEBUG("Try to acquire lock for pos " << u);
//			omp_set_lock(&mapLocks[u]);
//#endif
			double volN = 0.0;
#pragma omp atomic read
			volN = volNode[u];
//#ifdef _OPENMP
//			omp_unset_lock(&mapLocks[u]);
//			DEBUG("Released lock for pos " << u);
//#endif
			double delta = (omegaCut(u, D) - omegaCut(u, C)) / total + ((volClusterMinusNode(C, u) - volClusterMinusNode(D, u)) * volN) / (2 * total * total);
			return delta;
		};

	Luby luby; // independent set algorithm
	std::vector<bool> I;
	if (this->parallelism == "independent") {
		DEBUG("finding independent set");
		I = luby.run(G);
	}

#if 0
	std::vector<bool> indSet(n, true);
	bool change = false;

	int i = 0;
	do {
		change = false;

		G.parallelForNodes([&](node u) {
					index neighInI = none;
					G.forNeighborsOf(u, [&](node v) {
								if (indSet[v]) {
									if (neighInI != none) {
										indSet[v] = false;
										indSet[std::min(u, neighInI)] = false;
										neighInI = std::max(u, neighInI);
										change = true;
									}
									else {
										neighInI = v;
									}
								}
							});
				});
		++i;
		DEBUG("end of iteration " << i);

#ifdef DEBUG
		count summe = 0;
		G.forNodes([&](node u) {
					summe += indSet[u];
				});
		DEBUG("independent set size " << summe);
#endif
	} while (change);
#endif

	// begin pass
	int i = 0;
	bool change = false; // change in last iteration?
	do {
#pragma omp barrier
		i += 1;
		DEBUG("---> Louvain pass: iteration # " << i << ", parallelism: " << this->parallelism);
		change = false; // is clustering stable?

#pragma omp parallel for \
		shared(change, zeta, volCluster, volNode) \
		schedule(guided) // FIXME

		for (node u = 0; u < n; ++u) {
			// call here

		// try to improve modularity by moving a node to neighboring clusters
//		auto moveNode = [&](node u) {
//			std::cout << u << " ";

			cluster best = none;
			cluster C = none;
			cluster D = none;
			double deltaBest = -0.5;

//#ifdef _OPENMP
//			DEBUG("Try to acquire lock for pos " << u);
//			omp_set_lock(&mapLocks[u]);
//#endif
#pragma omp atomic read
			C = zeta[u];
//#ifdef _OPENMP
//			omp_unset_lock(&mapLocks[u]);
//			DEBUG("Released lock for pos " << u);
//#endif

//			TRACE("Processing neighborhood of node " << u << ", which is in cluster " << C);
			G.forNeighborsOf(u, [&](node v) {
//#ifdef _OPENMP
//				DEBUG("Try to acquire lock for pos " << v);
//				omp_set_lock(&mapLocks[v]);
//#endif
#pragma omp atomic read
				D = zeta[v];
//#ifdef _OPENMP
//				omp_unset_lock(&mapLocks[v]);
//				DEBUG("Released lock for pos " << v);
//#endif
//				TRACE("Neighbor " << v << ", which is still in cluster " << zeta[v]);
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
					DEBUG("Try to acquire lock for pos " << v);
					omp_set_lock(&mapLocks[v]);
#endif
					incidenceWeight[v][best] += w;
					incidenceWeight[v][C] -= w;
#ifdef _OPENMP
					omp_unset_lock(&mapLocks[v]);
					DEBUG("Released lock for pos " << v);
#endif
				});
//#ifdef _OPENMP
//				DEBUG("Try to acquire lock for pos " << u);
//				omp_set_lock(&mapLocks[u]);
//#endif
#pragma omp atomic write
				zeta[u] = best; // move to best cluster
#pragma omp atomic read
				volN = volNode[u];
//#ifdef _OPENMP
//				omp_unset_lock(&mapLocks[u]);
//				DEBUG("Released lock for pos " << u);
//#endif

#pragma omp atomic write
				change = true; // change to clustering has been made

#pragma omp atomic write
				this->anyChange = true; // indicate globally that clustering was modified

				// update the volume of the two clusters
//#ifdef _OPENMP
//				DEBUG("Try to acquire lock for pos " << C);
//				omp_set_lock(&mapLocks[C]);
//#endif
#pragma omp atomic update
				volCluster[C] -= volN;
//#ifdef _OPENMP
//				omp_unset_lock(&mapLocks[C]);
//				DEBUG("Released lock for pos " << C);
//				DEBUG("Try to acquire lock for pos " << best);
//				omp_set_lock(&mapLocks[best]);
//#endif
#pragma omp atomic update
				volCluster[best] += volN;
//#ifdef _OPENMP
//				omp_unset_lock(&mapLocks[best]);
//				DEBUG("Released lock for pos " << best);
//#endif
			}
//		};
		}

		// FIXME: uncomment
		// apply node movement according to parallelization strategy
//		if (this->parallelism == "none") {
//			G.forNodes(moveNode);
//		} else if (this->parallelism == "naive") {
//			G.parallelForNodes(moveNode);
//		} else if (this->parallelism == "naive-balanced") {
//			G.balancedParallelForNodes(moveNode);
//		} else if (this->parallelism == "independent") {
//			// try to move only the nodes in independent set
//			G.parallelForNodes([&](node u) {
//				if (I[u]) {
//					moveNode(u);
//				}
//			});
//		} else {
//			ERROR("unknown parallelization strategy: " << this->parallelism);
//			exit(1);
//		}

//		std::cout << std::endl;
#pragma omp barrier
	} while (change && i < MAX_LOUVAIN_ITERATIONS);

	// free lock mem
#pragma omp barrier
	DEBUG("about to destroy locks in Louvain pass");
#ifdef _OPENMP
	G.parallelForNodes([&](node u) {
				omp_destroy_lock(&mapLocks[u]);
			});
#endif

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

	int h = -1; // finest hierarchy level
	bool done = false; //

	Graph* graph = &G;

	do {
		h += 1; // begin new hierarchy level
		INFO("Louvain hierarchy level " << h);
		// one local optimization pass
		DEBUG("starting Louvain pass");
		Clustering clustering = this->pass(*graph); // FIXME
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

std::string Louvain::toString() const {
	std::stringstream strm;
	strm << "Louvain(" << "parallelism=" << this->parallelism << ")";
	return strm.str();
}

} /* namespace EnsembleClustering */
