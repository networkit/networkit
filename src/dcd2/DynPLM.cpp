/*
 * DynPLM.cpp
 *
 *  Created on: 03.01.2014
 *      Author: cls
 */

#include "DynPLM.h"

namespace	 NetworKit {

DynPLM::DynPLM(std::string prepStrategy, bool refine, double gamma,
		std::string par) :
		prepStrategy(prepStrategy), parallelism(par), refine(refine), gamma(gamma), plm2(refine, gamma, par) {

}

void DynPLM::update(std::vector<GraphEvent>& stream) {
	DEBUG("processing event stream");

	// turn a node into a singleton
	auto isolate = [&](node u) {
		TRACE("isolating " << u);
		zeta[u] = zeta.addCluster();
		TRACE("by creating singleton " << zeta[u]);
	};

	if (prepStrategy == "isolate") {
		for (GraphEvent ev : stream) {
			TRACE("event: " << ev.toString());
			switch (ev.type) {
				case GraphEvent::NODE_ADDITION : {
					zeta.append(ev.u);
					isolate(ev.u);
					break;
				}
				case GraphEvent::NODE_REMOVAL : {
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
	} else if (prepStrategy == "isolateNeighbors") {

		// turn a node into a singleton
		auto tryIsolate = [&](node u) {
			if (zeta.contains(u)) {
				zeta[u] = zeta.addCluster();
			}
		};

		for (GraphEvent ev : stream) {
			TRACE("event: " << ev.toString());
			TRACE("zeta: " << Aux::vectorToString(zeta.getVector()));
			switch (ev.type) {
				case GraphEvent::NODE_ADDITION : {
					// FIXME: segmentation fault 
					zeta.append(ev.u);
					isolate(ev.u);
					break;
				}
				case GraphEvent::NODE_REMOVAL : {
					zeta[ev.u] = none;
					break;
				}
				case GraphEvent::EDGE_ADDITION : {
					isolate(ev.u);
					isolate(ev.v);
					G->forNeighborsOf(ev.u, [&](node v) {
						tryIsolate(v);
					});
					G->forNeighborsOf(ev.v, [&](node v) {
						tryIsolate(v);
					});
					break;
				}
				case GraphEvent::EDGE_REMOVAL : {
					isolate(ev.u);
					isolate(ev.v);
					G->forNeighborsOf(ev.u, [&](node v) {
						tryIsolate(v);
					});
					G->forNeighborsOf(ev.v, [&](node v) {
						tryIsolate(v);
					});
					break;
				}
				case GraphEvent::EDGE_WEIGHT_UPDATE : {
					isolate(ev.u);
					isolate(ev.v);
					G->forNeighborsOf(ev.u, [&](node v) {
						tryIsolate(v);
					});
					G->forNeighborsOf(ev.v, [&](node v) {
						tryIsolate(v);
					});
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
	} else {
		ERROR("unknown prep strategy" << prepStrategy);
		throw std::runtime_error("unknown prep strategy");
	}
}

Clustering DynPLM::detect() {
	DEBUG("retrieving solution");
	return this->run(*this->G);
}

Clustering DynPLM::run(Graph& G) {
	INFO("calling run method on " << G.toString());

	edgeweight total = G.totalEdgeWeight();
	edgeweight divisor = (2 * total * total); // needed in modularity calculation


	// init community-dependent temporaries
	DEBUG("initializing community volume map");
	std::map<cluster, double> volCommunity; // a map to save memory
	zeta.forEntries([&](node u, cluster C) { 	// set volume for all communities
		volCommunity[C] += G.volume(u);
	});

	bool moved = false; // indicates whether any node has been moved


	// try to improve modularity by moving a node to neighboring clusters
	auto tryMove = [&](node u) {
		TRACE("trying to move node " << u);

		TRACE("zeta: " << Aux::vectorToString(zeta.getVector()));
		// collect edge weight to neighbor clusters
		std::map<cluster, edgeweight> affinity;
		G.forWeightedNeighborsOf(u, [&](node v, edgeweight weight) {
			if (u != v) {
				TRACE(u << " " << v);
				cluster C = zeta[v]; // FIXME: zeta[v] does not exist
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

		TRACE("Processing neighborhood of node " << u << ", which is in cluster " << C);
		G.forNeighborsOf(u, [&](node v) {
			D = zeta[v];
			if (D != C) { // consider only nodes in other clusters (and implicitly only nodes other than u)
				double delta = modGain(u, C, D);
				TRACE("mod gain: " << delta);
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
			// TRACE("node " << u << " moved");
			modUpdate(u, C, best);
			moved = true; // change to clustering has been made

		} else {
			// TRACE("node " << u << " not moved");
		}
	};

	// first move phase

	// apply node movement according to parallelization strategy
	DEBUG("starting move phase");
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
		std::pair<Graph, std::vector<node>> coarsened = plm2.coarsen(G, zeta);	// coarsen graph according to communitites
		Clustering zetaCoarse = plm2.run(coarsened.first);

		zeta = plm2.prolong(coarsened.first, zetaCoarse, G, coarsened.second); // unpack communities in coarse graph onto fine graph
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

