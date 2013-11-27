/*
 * MLPLM.cpp
 *
 *  Created on: 20.11.2013
 *      Author: cls
 */

#include "MLPLM.h"
#include <omp.h>
#include "../coarsening/ClusterContracter.h"
#include "../coarsening/ClusteringProjector.h"

#include <sstream>

namespace NetworKit {

MLPLM::MLPLM(std::string par, bool refine, double gamma) : parallelism(par), refine(refine), gamma(gamma){

}


Clustering MLPLM::run(Graph& G) {
	INFO("calling run method on " << G.toString());

	count z = G.upperNodeIdBound();


	// init communities to singletons
	Clustering zeta(z);
	zeta.allToSingletons();
	index o = zeta.upperBound();

	// init graph-dependent temporaries
	std::vector<double> volNode(z, 0.0);
	// $\omega(E)$
	edgeweight total = G.totalEdgeWeight();
	DEBUG("total edge weight: " << total);

	G.parallelForNodes([&](node u) { // calculate and store volume of each node
		volNode[u] += G.weightedDegree(u);
		volNode[u] += G.weight(u, u); // consider self-loop twice
		// TRACE("init volNode[" << u << "] to " << volNode[u]);
	});

	// init community-dependent temporaries
	std::vector<double> volCommunity(o, 0.0);
	zeta.parallelForEntries([&](node u, cluster C) { 	// set volume for all communities
		volCommunity[C] = volNode[u];
	});

	std::vector<std::map<cluster, edgeweight> > incidenceWeight(z); // For each node we store a map that maps from cluster ID to weight of edges to that cluster, this needs to be updated when a change occurs

	std::vector<omp_lock_t> mapLocks(z); // create and init locks
	G.parallelForNodes([&](node u) {
		omp_init_lock(&mapLocks[u]);
	});

	G.parallelForNodes([&](node u) {
		G.forWeightedEdgesOf(u, [&](node u, node v, edgeweight w) {
			cluster C = zeta[v];
			if (u != v) {
				incidenceWeight[u][C] += w;
			}
		});
	});
	bool moved = false; // indicates whether any node has been moved


	// try to improve modularity by moving a node to neighboring clusters
	auto tryMove = [&](node u) {
		TRACE("trying to move node " << u);
		// sub-functions

		// $\vol(C \ {x})$ - volume of cluster C excluding node x
		auto volCommunityMinusNode = [&](cluster C, node x) {
			double volC = 0.0;
			double volN = 0.0;
			volC = volCommunity[C];
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
			omp_set_lock(&mapLocks[u]);
			w = incidenceWeight[u][C];
			omp_unset_lock(&mapLocks[u]);
			return w;
		};

		auto modGain = [&](node u, cluster C, cluster D) {
			double volN = 0.0;
			volN = volNode[u];
			double delta = (omegaCut(u, D) - omegaCut(u, C)) / total + this->gamma * ((volCommunityMinusNode(C, u) - volCommunityMinusNode(D, u)) * volN) / (2 * total * total);
			TRACE("(" << omegaCut(u, D) << " - " << omegaCut(u, C) << ") / " << total << " + " << this->gamma << " * ((" << volCommunityMinusNode(C, u) << " - " << volCommunityMinusNode(D, u) << ") *" << volN << ") / 2 * " << (total * total));
			return delta;
		};


		auto modUpdate = [&](node u, cluster C, cluster D) {
			// update weight of edges to incident clusters
			G.forWeightedNeighborsOf(u, [&](node v, edgeweight w) {
				omp_set_lock(&mapLocks[v]);
				incidenceWeight[v][D] += w;
				incidenceWeight[v][C] -= w;
				omp_unset_lock(&mapLocks[v]);
			});

			double volN = 0.0;
			volN = volNode[u];
			// update the volume of the two clusters
			#pragma omp atomic update
			volCommunity[C] -= volN;
			#pragma omp atomic update
			volCommunity[D] += volN;
		};

		cluster best = none;
		cluster C = none;
		cluster D = none;
		double deltaBest = -1;

		C = zeta[u];

//			TRACE("Processing neighborhood of node " << u << ", which is in cluster " << C);
		G.forNeighborsOf(u, [&](node v) {
			D = zeta[v];
			if (D != C) { // consider only nodes in other clusters (and implicitly only nodes other than u)
				double delta = modGain(u, C, D);
				TRACE("mod gain: " << delta); // FIXME: all mod gains are negative
				if (delta > deltaBest) {
					deltaBest = delta;
					best = D;
				}
			}
		});

		TRACE("deltaBest=" << deltaBest); // FIXME: best mod gain is negative
		if (deltaBest > 0) { // if modularity improvement possible
			assert (best != C && best != none);// do not "move" to original cluster

			zeta[u] = best; // move to best cluster
			TRACE("node " << u << " moved");
			modUpdate(u, C, best);

			moved = true; // change to clustering has been made

		} else {
			TRACE("node " << u << " not moved");
		}
	};

	// first move phase

	// apply node movement according to parallelization strategy
	if (this->parallelism == "none") {
		G.forNodes(tryMove);
	} else if (this->parallelism == "simple") {
		G.parallelForNodes(tryMove);
	} else if (this->parallelism == "balanced") {
		G.balancedParallelForNodes(tryMove);
	} else {
		ERROR("unknown parallelization strategy: " << this->parallelism);
		throw std::runtime_error("unknown parallelization strategy");
	}


	if (moved) {
		std::pair<Graph, std::vector<node>> coarsened = coarsen(G, zeta);	// coarsen graph according to communitites
		Clustering zetaCoarse = run(coarsened.first);

		zeta = prolong(coarsened.first, zetaCoarse, G, coarsened.second); // unpack communities in coarse graph onto fine graph
		// refinement phase
		if (refine) {
			INFO("refinement phase");
			// reinit community-dependent temporaries
			o = zeta.upperBound();
			volCommunity.clear();
			volCommunity.resize(o, 0.0);
			zeta.parallelForEntries([&](node u, cluster C) { 	// set volume for all communities
				edgeweight volN = volNode[u];
				#pragma omp atomic update
				volCommunity[C] += volN;
			});

			incidenceWeight.clear(); // For each node we store a map that maps from cluster ID to weight of edges to that cluster, this needs to be updated when a change occurs
			incidenceWeight.resize(z);
			G.parallelForNodes([&](node u) {
				G.forWeightedEdgesOf(u, [&](node u, node v, edgeweight w) {
					cluster C = zeta[v];
					if (u != v) {
						#pragma omp atomic update
						incidenceWeight[u][C] += w;
					}
				});
			});

			// second move phase
			moved = false;
			// apply node movement according to parallelization strategy
			if (this->parallelism == "none") {
				G.forNodes(tryMove);
			} else if (this->parallelism == "simple") {
				G.parallelForNodes(tryMove);
			} else if (this->parallelism == "balanced") {
				G.balancedParallelForNodes(tryMove);
			} else {
				ERROR("unknown parallelization strategy: " << this->parallelism);
				throw std::runtime_error("unknown parallelization strategy");
			}
		}
	}

	return zeta;
}

std::string NetworKit::MLPLM::toString() const {
	std::string refined;
	if (refine) {
		refined = "refinement";
	} else {
		refined = "";
	}

	std::stringstream stream;
	stream << "MLPM(" << parallelism << "," << refined << ")";

	return stream.str();
}

std::pair<Graph, std::vector<node> > MLPLM::coarsen(const Graph& G, const Clustering& zeta) {
	Graph Gcoarse(0); // empty graph
	Gcoarse.markAsWeighted(); // Gcon will be a weighted graph

	std::vector<node> communityToMetaNode(zeta.upperBound(), none); // there is one meta-node for each community

	// populate map cluster -> meta-node
	G.forNodes([&](node v){
		cluster c = zeta[v];
		if (communityToMetaNode[c] == none) {
			communityToMetaNode[c] = Gcoarse.addNode(); // TODO: probably does not scale well, think about allocating ranges of nodes
		}
	});


	count z = G.upperNodeIdBound();
	std::vector<node> nodeToMetaNode(z);

	// set entries node -> supernode
	G.parallelForNodes([&](node v){
		nodeToMetaNode[v] = communityToMetaNode[zeta[v]];
	});


	// iterate over edges of G and create edges in Gcon or update edge and node weights in Gcon
	G.forWeightedEdges([&](node u, node v, edgeweight ew) {
		node su = nodeToMetaNode[u];
		node sv = nodeToMetaNode[v];
		if (zeta[u] == zeta[v]) {
			// add edge weight to supernode (self-loop) weight
			Gcoarse.setWeight(su, su, Gcoarse.weight(su, su) + ew);
		} else {
			// add edge weight to weight between two supernodes (or insert edge)
			Gcoarse.setWeight(su, sv, Gcoarse.weight(su, sv) + ew);
		}
	}); 

	return std::make_pair(Gcoarse, nodeToMetaNode);

}

Clustering MLPLM::prolong(const Graph& Gcoarse, const Clustering& zetaCoarse, const Graph& Gfine, std::vector<node> nodeToMetaNode) {
	Clustering zetaFine(Gfine.upperNodeIdBound());

	Gfine.forNodes([&](node v) {
		node mv = nodeToMetaNode[v];
		cluster cv = zetaCoarse[mv];
		zetaFine[v] = cv;
	});

	return zetaFine;
}

} /* namespace NetworKit */

