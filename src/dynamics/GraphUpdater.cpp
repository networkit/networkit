/*
 * GraphUpdater.cpp
 *
 *  Created on: 27.12.2013
 *      Author: cls
 */

#include "GraphUpdater.h"

namespace NetworKit {

GraphUpdater::GraphUpdater(Graph& G) : G(G) {
}

void GraphUpdater::update(std::vector<GraphEvent>& stream) {
	for (GraphEvent ev : stream) {
		switch (ev.type) {
			case GraphEvent::NODE_ADDITION : {
				node u = G.addNode();
				assert (u == ev.u);
				break;
			}
			case GraphEvent::NODE_REMOVAL : {
				G.removeNode(ev.u);
				break;
			}
			case GraphEvent::EDGE_ADDITION : {
				G.addEdge(ev.u, ev.v, ev.w);
				break;
			}
			case GraphEvent::EDGE_REMOVAL : {
				G.removeEdge(ev.u, ev.w);
				break;
			}
			case GraphEvent::EDGE_WEIGHT_UPDATE : {
				G.setWeight(ev.u, ev.v, ev.w);
				break;
			}
			case GraphEvent::TIME_STEP : {
				G.timeStep();
				break;
			}
			default: {
				throw std::runtime_error("unknown event type");
			}
		}
	}
}

} /* namespace NetworKit */



