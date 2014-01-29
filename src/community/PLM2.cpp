/*
 * MLPLM.cpp
 *
 *  Created on: 20.11.2013
 *      Author: cls
 */

#include "PLM2.h"
#include <omp.h>
#include "../coarsening/PartitionCoarsening.h"
 #include "../coarsening/ClusterContracter.h"
#include "../coarsening/ClusteringProjector.h"

#include <sstream>

namespace NetworKit {

PLM2::PLM2(bool refine, double gamma, std::string par, count maxIter) : parallelism(par), refine(refine), gamma(gamma), maxIter(maxIter) {

}


Clustering PLM2::run(Graph& G) {
	INFO("calling run method on " , G.toString());

	count z = G.upperNodeIdBound();


	// init communities to singletons
	Clustering zeta(z);
	zeta.allToSingletons();
	index o = zeta.upperBound();

	// init graph-dependent temporaries
	std::vector<double> volNode(z, 0.0);
	// $\omega(E)$
	edgeweight total = G.totalEdgeWeight();
	DEBUG("total edge weight: " , total);
	edgeweight divisor = (2 * total * total); // needed in modularity calculation

	G.parallelForNodes([&](node u) { // calculate and store volume of each node
		volNode[u] += G.weightedDegree(u);
		volNode[u] += G.weight(u, u); // consider self-loop twice
		// TRACE("init volNode[" , u , "] to " , volNode[u]);
	});

	// init community-dependent temporaries
	std::vector<double> volCommunity(o, 0.0);
	zeta.parallelForEntries([&](node u, cluster C) { 	// set volume for all communities
		volCommunity[C] = volNode[u];
	});

	// first move phase
	bool moved = false; // indicates whether any node has been moved in the last pass
	bool change = false; // indicates whether the communities have changed at all 

	// try to improve modularity by moving a node to neighboring clusters
	auto tryMove = [&](node u) {
		// TRACE("trying to move node " , u);

		// collect edge weight to neighbor clusters
		std::map<cluster, edgeweight> affinity;
		G.forWeightedNeighborsOf(u, [&](node v, edgeweight weight) {
			if (u != v) {
				cluster C = zeta[v];
				affinity[C] += weight;
			}
		});


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

		// // $\omega(u | C \ u)$
		// auto omegaCut = [&](node u, cluster C) {
		// 	return affinity[C];
		// };

		auto modGain = [&](node u, cluster C, cluster D) {
			double volN = 0.0;
			volN = volNode[u];
			double delta = (affinity[D] - affinity[C]) / total + this->gamma * ((volCommunityMinusNode(C, u) - volCommunityMinusNode(D, u)) * volN) / divisor; 
			//TRACE("(" , affinity[D] , " - " , affinity[C] , ") / " , total , " + " , this->gamma , " * ((" , volCommunityMinusNode(C, u) , " - " , volCommunityMinusNode(D, u) , ") *" , volN , ") / 2 * " , (total * total));
			return delta;
		};

		cluster best = none;
		cluster C = none;
		cluster D = none;
		double deltaBest = -1;

		C = zeta[u];

//			TRACE("Processing neighborhood of node " , u , ", which is in cluster " , C);
		G.forNeighborsOf(u, [&](node v) {
			D = zeta[v];
			if (D != C) { // consider only nodes in other clusters (and implicitly only nodes other than u)
				double delta = modGain(u, C, D);
				// TRACE("mod gain: " , delta); // FIXME: all mod gains are negative
				if (delta > deltaBest) {
					deltaBest = delta;
					best = D;
				}
			}
		});

		// TRACE("deltaBest=" , deltaBest); // FIXME: best mod gain is negative
		if (deltaBest > 0) { // if modularity improvement possible
			assert (best != C && best != none);// do not "move" to original cluster

			zeta[u] = best; // move to best cluster
			// TRACE("node " , u , " moved");

			// mod update
			double volN = 0.0;
			volN = volNode[u];
			// update the volume of the two clusters
			#pragma omp atomic update
			volCommunity[C] -= volN;
			#pragma omp atomic update
			volCommunity[best] += volN;

			moved = true; // change to clustering has been made

		} else {
			// TRACE("node " , u , " not moved");
		}
	};

	count iter = 0;
	do {
		iter += 1;
		moved = false;
		// apply node movement according to parallelization strategy
		if (this->parallelism == "none") {
			G.forNodes(tryMove);
		} else if (this->parallelism == "simple") {
			G.parallelForNodes(tryMove);
		} else if (this->parallelism == "balanced") {
			G.balancedParallelForNodes(tryMove);
		} else {
			ERROR("unknown parallelization strategy: " , this->parallelism);
			throw std::runtime_error("unknown parallelization strategy");
		}
		if (moved) change = true;
	} while (moved && (iter <= maxIter));
	INFO("iterations in move phase: ", iter);


	if (change) {
		DEBUG("nodes moved, so begin coarsening and recursive call");
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

			// second move phase
			iter = 0;
			do {
				iter += 1;
				moved = false;
				// apply node movement according to parallelization strategy
				if (this->parallelism == "none") {
					G.forNodes(tryMove);
				} else if (this->parallelism == "simple") {
					G.parallelForNodes(tryMove);
				} else if (this->parallelism == "balanced") {
					G.balancedParallelForNodes(tryMove);
				} else {
					ERROR("unknown parallelization strategy: " , this->parallelism);
					throw std::runtime_error("unknown parallelization strategy");
				}
			} while (moved && (iter <= maxIter));
			INFO("iterations in refinement move phase: ", iter);
		}
	}

	return zeta;
}

std::string NetworKit::PLM2::toString() const {
	std::string refined;
	if (refine) {
		refined = "refinement";
	} else {
		refined = "";
	}

	std::stringstream stream;
	stream << "PLM2(" << parallelism << "," << refined << ")";

	return stream.str();
}

std::pair<Graph, std::vector<node> > PLM2::coarsen(const Graph& G, const Clustering& zeta) {
	bool parallelCoarsening = false; // switch between parallel and sequential coarsening
	if (parallelCoarsening) {
		PartitionCoarsening parCoarsening;
		return parCoarsening.run(G, zeta);
	} else {
		ClusterContracter seqCoarsening;
		return seqCoarsening.run(G, zeta);
	}

}

Clustering PLM2::prolong(const Graph& Gcoarse, const Clustering& zetaCoarse, const Graph& Gfine, std::vector<node> nodeToMetaNode) {
	Clustering zetaFine(Gfine.upperNodeIdBound());

	Gfine.forNodes([&](node v) {
		node mv = nodeToMetaNode[v];
		cluster cv = zetaCoarse[mv];
		zetaFine[v] = cv;
	});

	return zetaFine;
}

} /* namespace NetworKit */

