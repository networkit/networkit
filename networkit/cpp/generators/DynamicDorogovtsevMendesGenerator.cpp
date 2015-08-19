/*
 * DynamicDorogovtsevMendesGenerator.cpp
 *
 *  Created on: 03.02.2014
 *      Author: cls
 */

#include "DynamicDorogovtsevMendesGenerator.h"

namespace NetworKit {

DynamicDorogovtsevMendesGenerator::DynamicDorogovtsevMendesGenerator() : initial(true), u(0) {

}

std::vector<GraphEvent> DynamicDorogovtsevMendesGenerator::generate(count nSteps) {
	
	std::vector<GraphEvent> stream;

	if (initial) {
		node s1 = u++;
		stream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, s1));
		node s2 = u++;
		stream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, s2));
		node s3 = u;
		stream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, s3));
		edges.push_back({s1, s2});
		stream.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, s1, s2));
		edges.push_back({s2, s3});
		stream.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, s2, s3));
		edges.push_back({s3, s1});
		stream.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, s3, s1));

		stream.push_back(GraphEvent(GraphEvent::TIME_STEP));
		initial = false;
	}

	for (index i = 0; i < nSteps; ++i) {
		++u; // new node
		stream.push_back(GraphEvent(GraphEvent::NODE_ADDITION, u));
		// select random edge
		index e = Aux::Random::integer(edges.size() - 1);
		node s = edges[e].first;
		node t = edges[e].second;
		edges.push_back({s,u});
		edges.push_back({t,u});
		// connect node
		stream.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, u, s));
		stream.push_back(GraphEvent(GraphEvent::EDGE_ADDITION, u, t));

		stream.push_back(GraphEvent(GraphEvent::TIME_STEP));
	}

	return stream;
}


} /* namespace NetworKit */


