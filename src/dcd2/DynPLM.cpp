/*
 * DynPLM.cpp
 *
 *  Created on: 03.01.2014
 *      Author: cls
 */

#include "DynPLM.h"

namespace	 NetworKit {

DynPLM::DynPLM(Graph& G, bool refine, double gamma, std::string par) : DynCommunityDetector(G), parallelism(par), refine(refine), gamma(gamma) {
	if (G.numberOfNodes() != 0) {
		throw std::runtime_error("DynPLM must be initialized with an empty graph");
	}


}

void DynPLM::process(std::vector<GraphEvent>& stream) {
	DEBUG("processing event stream");

	auto isolate = [&](node u) {
		edgeweight vol = G.volume(u);
		// volCommunity[zeta[u]] -= vol;
		zeta[u] = zeta.addCluster();
		// volCommunity[zeta[u]] = vol;
	};

	for (GraphEvent ev : stream) {
		TRACE("event: " << ev.toString());
		switch (ev.type) {
			case GraphEvent::NODE_ADDITION : {
				zeta.append(ev.u);
				zeta[ev.u] = zeta.addCluster();
				// volCommunity[zeta[ev.u]] = G.volume(ev.u);
				break;
			}
			case GraphEvent::NODE_REMOVAL : {
				// volCommunity[zeta[ev.u]] -= G.volume(ev.u);
				zeta[ev.u] = none;
				break;
			}
			case GraphEvent::EDGE_ADDITION : {
				isolate(ev.u);
				isolate(ev.v);
				break;
			}
			case GraphEvent::EDGE_REMOVAL : {
				isolate(ev.u);
				isolate(ev.v);
				break;
			}
			case GraphEvent::EDGE_WEIGHT_UPDATE : {
				isolate(ev.u);
				isolate(ev.v);
				break;
			}
			case GraphEvent::TIME_STEP : {
				break;
			}
			default: {
				throw std::runtime_error("unknown event type");
			}
		}
	} // end event loop
}

Clustering DynPLM::retrieve() {
	DEBUG("retrieving solution");
	return this->run(this->G);
}

Clustering DynPLM::run(Graph& G) {
	INFO("calling run method on " << G.toString());

	edgeweight total = G.totalEdgeWeight();
	edgeweight divisor = (2 * total * total); // needed in modularity calculation



	// init community-dependent temporaries
	std::map<cluster, double> volCommunity;
	// std::vector<double> volCommunity(zeta.upperBound(), 0.0);
	zeta.forEntries([&](node u, cluster C) { 	// set volume for all communities
		volCommunity[C] += G.volume(u);
	});

	bool moved = false; // indicates whether any node has been moved


	// try to improve modularity by moving a node to neighboring clusters
	auto tryMove = [&](node u) {
		// TRACE("trying to move node " << u);

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
			volC = volCommunity[C];
			if (zeta[x] == C) {
				return volC - G.volume(x);
			} else {
				return volC;
			}
		};

		auto modGain = [&](node u, cluster C, cluster D) {
			double delta = (affinity[D] - affinity[C]) / total + this->gamma * ((volCommunityMinusNode(C, u) - volCommunityMinusNode(D, u)) * G.volume(u)) / divisor;
			//TRACE("(" << affinity[D] << " - " << affinity[C] << ") / " << total << " + " << this->gamma << " * ((" << volCommunityMinusNode(C, u) << " - " << volCommunityMinusNode(D, u) << ") *" << volN << ") / 2 * " << (total * total));
			return delta;
		};


		auto modUpdate = [&](node u, cluster C, cluster D) {
			double volN = 0.0;
			volN = G.volume(u);
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
				// TRACE("mod gain: " << delta); // FIXME: all mod gains are negative
				if (delta > deltaBest) {
					deltaBest = delta;
					best = D;
				}
			}
		});

		// TRACE("deltaBest=" << deltaBest); // FIXME: best mod gain is negative
		if (deltaBest > 0) { // if modularity improvement possible
			assert (best != C && best != none);// do not "move" to original cluster

			zeta[u] = best; // move to best cluster
			// TRACE("node " << u << " moved");
			modUpdate(u, C, best);

			moved = true; // change to clustering has been made

		} else {
			// TRACE("node " << u << " not moved");
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
		DEBUG("nodes moved, so begin coarsening and recursive call");
		std::pair<Graph, std::vector<node>> coarsened = PLM2::coarsen(G, zeta);	// coarsen graph according to communitites
		Clustering zetaCoarse = run(coarsened.first);

		zeta = PLM2::prolong(coarsened.first, zetaCoarse, G, coarsened.second); // unpack communities in coarse graph onto fine graph
		// refinement phase
		if (refine) {
			INFO("refinement phase");
			// reinit community-dependent temporaries
			cluster o = zeta.upperBound();
			volCommunity.clear();
			zeta.forEntries([&](node u, cluster C) { 	// set volume for all communities
				volCommunity[C] += G.volume(u);
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


} /* namespace NetworKit */

