/*
 * DynPLP.cpp
 *
 *  Created on: 03.01.2014
 *      Author: cls
 */

#include "DynPLP.h"

namespace NetworKit {

DynPLP::DynPLP(std::string prepStrategy, count theta) : prepStrategy(prepStrategy), updateThreshold(theta) {

}

void DynPLP::update(std::vector<GraphEvent>& stream) {
	DEBUG("processing event stream of length " , stream.size());
	auto isolate = [&](node u) {
		zeta.toSingleton(u);
		activeNodes[u] = true;
	};

	if (prepStrategy == "isolate") {
		for (GraphEvent ev : stream) {
			// TODO: remove trace
			// TRACE("event: " , ev.toString());
			switch (ev.type) {
				case GraphEvent::NODE_ADDITION : {
					index u = zeta.extend(); //actually returns an index
					zeta.toSingleton(u); //ev.u
					activeNodes.push_back(true);
					break;
				}
				case GraphEvent::NODE_REMOVAL : {
					zeta.remove(ev.u);
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

		auto tryIsolate = [&](node u) {
			if (zeta.contains(u)) {
				// because the graph can be in the future, data structures do not necessarily know neighbor node u
				zeta.toSingleton(u);
				activeNodes[u] = true;
			}
		};

		for (GraphEvent ev : stream) {
			// TRACE("event: " , ev.toString());
			switch (ev.type) {
				case GraphEvent::NODE_ADDITION : {
					index u = zeta.extend(); 
					activeNodes.push_back(true);
					isolate(ev.u); //u
					break;
				}
				case GraphEvent::NODE_REMOVAL : {
					zeta.remove(ev.u);
					break;
				}
				case GraphEvent::EDGE_ADDITION : {
					isolate(ev.u);
					isolate(ev.v);
					G->forNeighborsOf(ev.u, [&](node v) {
						// TRACE("neighbor " , v);
						tryIsolate(v);
					});
					G->forNeighborsOf(ev.v, [&](node v) {
						// TRACE("neighbor" , v);
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
		ERROR("unknown prep strategy: " , prepStrategy);
		throw std::runtime_error("unknown prep strategy");
	}
}

Partition DynPLP::detect(bool restart) {
	DEBUG("retrieving communities");

	if (restart) {
		zeta.allToSingletons();
	}

	nIterations = 0; // number of iterations

	count nUpdated;

	// propagate labels
	do {
		nUpdated = 0; // number of nodes which have been updated in last iteration
		nIterations++;
		// TRACE("iteration " , nIterations);

		auto propagate = [&](node v) {
			// TRACE("scanning node " , v);
			if ((activeNodes[v]) && (G->degree(v) > 0)) {
				// TRACE("processing node " , v);
				std::map<label, double> labelWeights; // neighborLabelCounts maps label -> frequency in the neighbors
				// weigh the labels in the neighborhood of v
				G->forWeightedNeighborsOf(v, [&](node w, edgeweight weight) {
					labelWeights[zeta[w]] += weight; // add weight of edge {v, w}
				});

				// get dominant label
				label dominant = std::max_element(labelWeights.begin(), labelWeights.end(), [](const std::pair<label, edgeweight>& p1, const std::pair<label, edgeweight>& p2) {return p1.second < p2.second;})->first;

				if (zeta[v] != dominant) { // UPDATE
					zeta[v] = dominant;
					nUpdated += 1; // TODO: atomic update?
					G->forNeighborsOf(v, [&](node u) {
						activeNodes[u] = true;
					});
				} else {
					activeNodes[v] = false;
				}
			}
		};

		G->balancedParallelForNodes(propagate);

		DEBUG("nodes updated: " , nUpdated);

	} while (nUpdated > this->updateThreshold);

	return zeta;

}

} /* namespace NetworKit */
