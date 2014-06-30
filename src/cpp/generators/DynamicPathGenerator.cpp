/*
 * DynamicPathGenerator.cpp
 *
 *  Created on: 14.01.2014
 *      Author: cls
 */

#include "DynamicPathGenerator.h"

namespace NetworKit {

std::vector<GraphEvent> DynamicPathGenerator::generate(count nSteps) {

	std::vector<GraphEvent> stream;
	node u = G.addNode();
	stream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, u));

	count step = 0;
	while (step < nSteps) {
		node v = G.addNode();
		stream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, v));
		G.addEdge(u, v);
		stream.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, u, v, 1.0));
		u = v;
		stream.push_back(GraphEvent(GraphEvent::TIME_STEP));
		step += 1;
	}
	return stream;
}

} /* namespace NetworKit */
