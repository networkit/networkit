/*
 * IsolateAffectedNodes.cpp
 *
 *  Created on: 27.12.2013
 *      Author: cls
 */

#include "IsolateAffectedNodes.h"

namespace NetworKit {

void IsolateAffectedNodes::process(Partition& communities, std::vector<GraphEvent> eventStream) {

	for (GraphEvent ev : eventStream) {
		switch (ev.type) {
			case GraphEvent::NODE_ADDITION : {
				node u = communities.extend();
				assert (u == ev.u);
				communities.toSingleton(ev.u);
				break;
			}
			case GraphEvent::NODE_REMOVAL : {
				assert (communities.contains(ev.u));
				communities.remove(ev.u);
				break;
			}
			case GraphEvent::EDGE_ADDITION : {
				assert (communities.contains(ev.u));
				assert (communities.contains(ev.v));
				communities.toSingleton(ev.u);
				communities.toSingleton(ev.v);
				break;
			}
			case GraphEvent::EDGE_REMOVAL : {
				assert (communities.contains(ev.u));
				assert (communities.contains(ev.v));
				communities.toSingleton(ev.u);
				communities.toSingleton(ev.v);
				break;
			}
			case GraphEvent::EDGE_WEIGHT_UPDATE : {
				assert (communities.contains(ev.u));
				assert (communities.contains(ev.v));
				communities.toSingleton(ev.u);
				communities.toSingleton(ev.v);
				break;
			}
			case GraphEvent::TIME_STEP : {
				// time step is ignored
				break;
			}
			default: {
				throw std::runtime_error("unknown event type");
			}
		}
	}

} /* namespace NetworKit */


}
